"""
02_mechanical_properties.py
===========================
Computes the tensile properties of every film from its stress-strain curve and
joins them with the film thickness.

Inputs
    raw_data/tensile_curves/tensile_curves.xlsx   (sheet "curves", one row per point)
    processed_data/thickness/thickness_films.csv    (from 01_thickness.py)
Outputs (processed_data/mech_properties/)
    mechanical_properties_films.csv      one row per film (108 films)
    mechanical_properties_architectures.csv   mean, SD and CV per architecture (36)
    mechprop_paper_format.xlsx           same film table with the column names used in
                                         the original analysis notebooks (mechprop.xlsx)
Figures (figures/)
    stress_strain_{family}.png           all curves of each family, one panel per film number

Definitions (per curve; strain in %, stress in MPa)
    E_slope_MPa_per_pct  slope of a linear fit to the initial part of the curve
                         (stress range in E_FIT_RANGE, below)
    E_MPa                E_slope * 100   (true Young's modulus in MPa)
    UTS_MPa              maximum stress
    strain_at_break_pct  maximum strain
    toughness_MPa_pct    trapezoidal area under the curve (paper convention)
    toughness_MJ_m3      toughness_MPa_pct / 100
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import lbl_common as C

# Young's modulus: stress range (MPa) of the linear fit in the initial part of the
# stress-strain curve, for each set of films
E_FIT_RANGE = {"P": (0.0, 5.0), "S": (0.0, 5.0), "NS": (2.5, 6.0), "SS": (2.5, 6.0)}

all_curves = C.read_tensile_curves()                     # {(architecture, film): (strain, stress)}
rows, curves_all = [], {}
for fam in C.FAMILIES:
    curves = {k: v for k, v in all_curves.items() if k[0].rstrip("0123456789") == fam}
    curves_all[fam] = curves
    for (arch, film), (strain, stress) in curves.items():
        props = C.curve_properties(strain, stress, E_FIT_RANGE[fam])
        rows.append({"family": fam, "architecture": arch, "film": film, **props})

films = pd.DataFrame(rows)
thick = pd.read_csv(C.PROC_THICK / "thickness_films.csv")[["architecture", "film", "thickness_um"]]
films = films.merge(thick, on=["architecture", "film"], how="left", validate="one_to_one")
films["E_MPa"] = films.E_slope_MPa_per_pct * 100
films["toughness_MJ_m3"] = films.toughness_MPa_pct / 100
films = films.sort_values(["architecture", "film"],
                          key=lambda s: s.map(C.arch_sort_key) if s.name == "architecture" else s)
cols = ["family", "architecture", "film", "thickness_um", "E_slope_MPa_per_pct", "E_MPa",
        "UTS_MPa", "strain_at_break_pct", "toughness_MPa_pct", "toughness_MJ_m3"]
films = films[cols].reset_index(drop=True)
assert films.thickness_um.notna().all(), "a film has no thickness"
films.to_csv(C.PROC_MECH / "mechanical_properties_films.csv", index=False)

# ---------------------------------------------------------------- per architecture
props = ["thickness_um", "E_slope_MPa_per_pct", "E_MPa", "UTS_MPa", "strain_at_break_pct",
         "toughness_MPa_pct", "toughness_MJ_m3"]
g = films.groupby(["family", "architecture"], sort=False)[props]
arch = g.mean().add_suffix("_mean").join(g.std(ddof=1).add_suffix("_sd"))
for p in props:
    arch[p + "_cv_pct"] = 100 * arch[p + "_sd"] / arch[p + "_mean"]
arch.insert(0, "n_films", g.size())
arch = arch[["n_films"] + [f"{p}_{s}" for p in props for s in ("mean", "sd", "cv_pct")]].reset_index()
arch.to_csv(C.PROC_MECH / "mechanical_properties_architectures.csv", index=False)

# ---------------------------------------------------------------- legacy format
legacy = films.rename(columns={"architecture": "Sample", "family": "Type", "film": "Rep",
                               "thickness_um": "Thickness_mm", "E_slope_MPa_per_pct": "Young_Modulus_MPa",
                               "strain_at_break_pct": "Elongation_%", "toughness_MPa_pct": "Toughness"})
legacy = legacy[["Sample", "Type", "Rep", "Thickness_mm", "Young_Modulus_MPa", "UTS_MPa", "Elongation_%", "Toughness"]]
legacy.to_excel(C.PROC_MECH / "mechprop_paper_format.xlsx", index=False)

# ---------------------------------------------------------------- figures
for fam, curves in curves_all.items():
    archs = sorted({a for a, _ in curves}, key=C.arch_sort_key)
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(archs), 2)))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for color, a in zip(cmap, archs):
        for film in (1, 2, 3):
            if (a, film) in curves:
                x, y = curves[(a, film)]
                axes[film - 1].plot(x, y, color=color, lw=1.2, label=a)
    for k, ax in enumerate(axes):
        ax.set_title(f"Film {k + 1}")
        ax.set_xlabel("Strain (%)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Stress (MPa)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=min(len(archs), 8), frameon=False, fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    fig.savefig(C.FIG / f"stress_strain_{fam}.png", dpi=300)
    plt.close(fig)

print(films.groupby("family").size().rename("films").to_string())
print(f"Written mechanical_properties_films.csv, "
      f"mechanical_properties_architectures.csv, mechprop_paper_format.xlsx")
