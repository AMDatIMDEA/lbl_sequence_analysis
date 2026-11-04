"""
04_statistics.py
================
Group statistics for the S, SS and NS libraries (Supporting Information, Sections 6-8).

Input   processed_data/mechanical_properties_films.csv   (from 02)
        processed_data/mechanical_properties_architectures.csv
Outputs (processed_data/statistics/)
    descriptive_by_family.csv     mean, SD, n and CV of every property per family (film level)
    welch_anova.csv               Welch one-way ANOVA S vs SS vs NS (film level), plus Shapiro and Levene
    tukey_<property>.csv          Tukey HSD pairwise comparisons (film level)
    pearson_architecture_means.csv  Pearson correlation matrix of architecture means (S+SS+NS)
Figures (figures/)
    FigS3_mechanical_properties.png   films (small) and architecture means (large) per family
    FigS4_pearson_matrix.png          correlation matrix of the architecture means

Properties are analysed in the units used in the paper (E as slope in MPa/%, toughness in
MPa*%). Correlations and p values do not depend on that choice of units.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import shapiro, levene, tukey_hsd, f as f_dist

import lbl_common as C

OUT = C.PROC / "statistics"
OUT.mkdir(exist_ok=True)
films = pd.read_csv(C.PROC / f"mechanical_properties_films.csv")
arch = pd.read_csv(C.PROC / f"mechanical_properties_architectures.csv")

ORDER = ["S", "SS", "NS"]
PROPS = {"thickness_um": "Thickness (µm)", "E_slope_MPa_per_pct": "Young's modulus (MPa/%)",
         "UTS_MPa": "UTS (MPa)", "strain_at_break_pct": "Elongation at break (%)",
         "toughness_MPa_pct": "Toughness (MPa·%)"}
d = films[films.family.isin(ORDER)]


def welch_anova(groups):
    """Welch's one-way ANOVA for unequal variances (same result as
    statsmodels.stats.oneway.anova_oneway(..., use_var="unequal"), used originally)."""
    k = len(groups)
    n = np.array([len(g) for g in groups], float)
    m = np.array([np.mean(g) for g in groups])
    v = np.array([np.var(g, ddof=1) for g in groups])
    w = n / v
    mw = (w * m).sum() / w.sum()
    a = (w * (m - mw) ** 2).sum() / (k - 1)
    lam = ((1 - w / w.sum()) ** 2 / (n - 1)).sum()
    b = 1 + 2 * (k - 2) / (k ** 2 - 1) * lam
    F = a / b
    df2 = (k ** 2 - 1) / (3 * lam)
    return F, f_dist.sf(F, k - 1, df2)

# ---------------------------------------------------------------- descriptive
desc = (d.groupby("family")[list(PROPS)].agg(["mean", "std", "count"]).stack(level=0, future_stack=True)
        .rename_axis(["family", "property"]).reset_index())
desc["cv_pct"] = 100 * desc["std"] / desc["mean"]
desc.to_csv(OUT / f"descriptive_by_family.csv", index=False)

# ---------------------------------------------------------------- tests
rows = []
for p in PROPS:
    groups = [d.loc[d.family == g, p].dropna() for g in ORDER]
    F, pw = welch_anova(groups)
    row = {"property": p, "welch_F": F, "welch_p": pw, "levene_p": levene(*groups).pvalue}
    row.update({f"shapiro_p_{g}": shapiro(x).pvalue for g, x in zip(ORDER, groups)})
    rows.append(row)
    t = tukey_hsd(*groups)                      # Tukey HSD (scipy); same as statsmodels pairwise_tukeyhsd
    ci = t.confidence_interval(0.95)
    tk = [{"group1": ORDER[i], "group2": ORDER[j], "meandiff": groups[j].mean() - groups[i].mean(),
           "p_adj": t.pvalue[i, j], "lower": -ci.high[i, j], "upper": -ci.low[i, j], "reject": t.pvalue[i, j] < 0.05}
          for i in range(3) for j in range(i + 1, 3)]
    pd.DataFrame(tk).to_csv(OUT / f"tukey_{p}.csv", index=False)
pd.DataFrame(rows).to_csv(OUT / f"welch_anova.csv", index=False)

# ---------------------------------------------------------------- Pearson of means
A = arch[arch.family.isin(ORDER)].set_index("architecture")[[p + "_mean" for p in PROPS]]
A.columns = list(PROPS.values())
corr = A.corr(method="pearson")
corr.to_csv(OUT / f"pearson_architecture_means.csv")

# ---------------------------------------------------------------- figures
colors = {"S": "#2a78d6", "SS": "#c9551f", "NS": "#3a9a5b"}
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, p in zip(axes.ravel(), ["E_slope_MPa_per_pct", "UTS_MPa", "strain_at_break_pct", "toughness_MPa_pct"]):
    x0 = 0
    for g in ORDER:
        names = sorted(d.loc[d.family == g, "architecture"].unique(), key=C.arch_sort_key)
        for k, a in enumerate(names):
            v = d.loc[d.architecture == a, p]
            ax.scatter(np.full(len(v), x0 + k), v, s=12, color=colors[g], alpha=0.5)
            ax.scatter(x0 + k, v.mean(), s=45, color=colors[g], edgecolor="black", lw=0.5, zorder=3)
        ax.text(x0 + (len(names) - 1) / 2, 1.02, g, transform=ax.get_xaxis_transform(), ha="center")
        x0 += len(names) + 1
    ax.set_ylabel(PROPS[p]); ax.set_xticks([]); ax.grid(alpha=0.2)
fig.tight_layout(); fig.savefig(C.FIG / f"FigS3_mechanical_properties.png", dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.5, 5.5))
im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(corr)), corr.columns, rotation=40, ha="right"); ax.set_yticks(range(len(corr)), corr.index)
for i in range(len(corr)):
    for j in range(len(corr)):
        ax.text(j, i, f"{corr.iat[i, j]:.2f}", ha="center", va="center", fontsize=9)
fig.colorbar(im, ax=ax, label="Pearson r"); fig.tight_layout()
fig.savefig(C.FIG / f"FigS4_pearson_matrix.png", dpi=300); plt.close(fig)

print(pd.DataFrame(rows)[["property", "welch_F", "welch_p"]].round(4).to_string(index=False))
