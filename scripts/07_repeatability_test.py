"""
07_repeatability_test.py
========================
Repeatability of the automated assembly: 10 independently prepared Alg/Chi (AX) films.

Inputs  raw_data/repeatability_test/repeatability_tensile_curves.xlsx, sheet "curves"
            one row per recorded point: Film, Point, Time (min), Temperature (°C), Strain (%), Stress (MPa)
        raw_data/repeatability_test/repeatability_thickness.xlsx, sheet "thickness_readings"
            one row per reading: Film, Reading, Thickness (mm), Thickness (µm); 8 readings per film
Outputs (processed_data/repeatability_test/)
    repeatability_films.csv     properties of every film
    repeatability_summary.csv   mean, SD and CV across films, with and without film 7
Figures figures/repeatability_curves.png

Two Young's-modulus definitions are reported:
    E_10_40UTS_MPa   slope between 10 % and 40 % of the UTS on the loading branch, x100
                     (definition used in the original repeatability analysis, rep_mechprop.xlsx)
    E_slope_0_5MPa   slope in the 0-5 MPa stress range, MPa/% (same fit as for the P films of the
                     library, which include Alg/Chi), for direct comparison with the main data set
Film 7 was treated as an outlier in the original analysis. Note that one of its thickness
readings is 0.150 mm (ten times the others), which looks like a transcription error.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lbl_common as C

OUT = C.PROC_REP

curves = pd.read_excel(C.RAW / "repeatability_test" / "repeatability_tensile_curves.xlsx", sheet_name="curves")
th = pd.read_excel(C.RAW / "repeatability_test" / "repeatability_thickness.xlsx", sheet_name="thickness_readings")

rows = []
fig, ax = plt.subplots(figsize=(6, 4.5))
for k in range(1, 11):
    c = curves[curves.Film == k].sort_values("Point")
    c = c[c["Strain (%)"].notna() & c["Stress (MPa)"].notna()]
    x, y = c["Strain (%)"].to_numpy(), c["Stress (MPa)"].to_numpy()
    uts, imax = y.max(), y.argmax()
    m = (y >= 0.1 * uts) & (y <= 0.4 * uts) & (np.arange(len(y)) <= imax)
    t = th.loc[th.Film == k, "Thickness (mm)"].dropna().to_numpy() * 1000      # mm -> um
    rows.append({"film": k,
                 "thickness_um": t.mean(), "thickness_max_um": t.max(), "n_thickness_readings": len(t),
                 "E_10_40UTS_MPa": np.polyfit(x[m], y[m], 1)[0] * 100,
                 "E_slope_0_5MPa_MPa_per_pct": C.young_slope(x, y, 0, 5),
                 "UTS_MPa": uts, "strain_at_break_pct": x.max(),
                 "toughness_MPa_pct": np.trapezoid(y, x), "toughness_MJ_m3": np.trapezoid(y, x) / 100})
    ax.plot(x, y, lw=1.2, label=f"Film {k}", ls="--" if k == 7 else "-")
ax.set_xlabel("Strain (%)"); ax.set_ylabel("Stress (MPa)"); ax.grid(alpha=0.3)
ax.legend(fontsize=7, ncol=2, frameon=False); fig.tight_layout()
fig.savefig(C.FIG / "repeatability_curves.png", dpi=300); plt.close(fig)

films = pd.DataFrame(rows)
films.to_csv(OUT / "repeatability_films.csv", index=False)


def summary(df, label):
    s = df.drop(columns="film")
    return pd.DataFrame({"set": label, "property": s.columns, "mean": s.mean().values,
                         "sd": s.std(ddof=1).values, "cv_pct": (100 * s.std(ddof=1) / s.mean()).values})


summ = pd.concat([summary(films, "all 10 films"), summary(films[films.film != 7], "without film 7")])
summ.to_csv(OUT / "repeatability_summary.csv", index=False)
print(summ.round(3).to_string(index=False))
