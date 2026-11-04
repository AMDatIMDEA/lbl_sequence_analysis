"""
05_parent_variant.py
====================
Compares every NS and SS variant with its parent S architecture (S1-S8)
(Supporting Information, Section 8; main text, sequence-perturbation analysis).

Inputs  processed_data/mechanical_properties_architectures.csv   (from 02)
        processed_data/parent_variant_ncd.csv                    (from 03)
Outputs (processed_data/statistics/)
    parent_variant_changes.csv   signed and relative change of every property, with NCD
    SS_wilcoxon.csv              paired Wilcoxon signed-rank test, S1-S8 vs SS1-SS8
    NS_spearman_vs_NCD.csv       Spearman correlation between NCD and |relative change| (S->NS)

Changes use architecture means:  signed % = 100 (variant - parent) / |parent|
"""
import pandas as pd
from scipy.stats import wilcoxon, spearmanr
import lbl_common as C

OUT = C.PROC / "statistics"
OUT.mkdir(exist_ok=True)
arch = pd.read_csv(C.PROC / f"mechanical_properties_architectures.csv").set_index("architecture")
ncd = pd.read_csv(C.PROC / "parent_variant_ncd.csv")
PROPS = ["thickness_um", "E_slope_MPa_per_pct", "E_MPa", "UTS_MPa", "strain_at_break_pct",
         "toughness_MPa_pct", "toughness_MJ_m3"]

rows = []
for _, r in ncd.iterrows():
    row = {"parent": r.parent, "variant": r.variant, "type": r.type, "NCD": r.NCD,
           "swapped_positions": r.swapped_positions}
    for p in PROPS:
        a, b = arch.loc[r.parent, p + "_mean"], arch.loc[r.variant, p + "_mean"]
        row[f"{p}_parent"] = a
        row[f"{p}_variant"] = b
        row[f"{p}_delta"] = b - a
        row[f"{p}_change_pct"] = 100 * (b - a) / abs(a)
    rows.append(row)
ch = pd.DataFrame(rows)
ch.to_csv(OUT / f"parent_variant_changes.csv", index=False)

main = ["thickness_um", "E_slope_MPa_per_pct", "UTS_MPa", "strain_at_break_pct", "toughness_MPa_pct"]
ss = ch[ch.type == "SS"]
pd.DataFrame([{"property": p, "wilcoxon_p": wilcoxon(ss[p + "_parent"], ss[p + "_variant"]).pvalue,
               "median_abs_change_pct": ss[p + "_change_pct"].abs().median()} for p in main]) \
  .to_csv(OUT / f"SS_wilcoxon.csv", index=False)

ns = ch[ch.type == "NS"]
sp = []
for p in main:
    rho, pv = spearmanr(ns.NCD, ns[p + "_change_pct"].abs())
    sp.append({"property": p, "spearman_rho": rho, "p": pv})
pd.DataFrame(sp).to_csv(OUT / f"NS_spearman_vs_NCD.csv", index=False)

show = ch[["parent", "variant", "NCD", "E_slope_MPa_per_pct_change_pct", "UTS_MPa_change_pct", "toughness_MPa_pct_change_pct"]]
print(show.round(2).to_string(index=False))
