"""
03_sequence_analysis.py
=======================
Sequence descriptors of the 36 architectures used in the paper.

Input
    raw_data/sequences/sequences.txt     deposition order of every architecture
                                         (A=ALG, B=PEC, X=CHI, Y=GEL; first element = first layer)
Outputs (processed_data/)
    sequence_composition.csv        number of layers and of each polymer per architecture
    ncd_matrix.csv                  pairwise Normalized Compression Distance (zlib), Eq. (1) of the SI
    parent_variant_ncd.csv          NCD and swapped positions for S->NS and S->SS pairs
    transition_counts.csv           the 8 directed nearest-neighbour transitions per architecture
    transition_fractions.csv        same, divided by the number of transitions (n_layers - 1)
    pairing_bias.csv                B_A->C and B_C->A descriptors (SI Section 9.2)
    parent_variant_dpair.csv        D_pair between each S parent and its NS/SS variant (SI Section 9.3)

NCD(x, y) = [C(xy) - min{C(x), C(y)}] / max{C(x), C(y)}
where C(.) is the length in bytes of the zlib-compressed sequence string
(default compression level), as in the original library-design notebook.
"""
import zlib
import numpy as np
import pandas as pd
import lbl_common as C

seqs = C.read_sequences()
archs = sorted([a for a in seqs if a.rstrip("0123456789") in C.FAMILIES], key=C.arch_sort_key)
text = {a: "".join(seqs[a]) for a in archs}

# ---------------------------------------------------------------- composition
comp = pd.DataFrame([{"architecture": a, "family": a.rstrip("0123456789"), "n_layers": len(seqs[a]),
                      **{NAME: seqs[a].count(code) for code, NAME in C.NAME.items()}} for a in archs])
comp.to_csv(C.PROC / "sequence_composition.csv", index=False)

# ---------------------------------------------------------------- NCD
def csize(s):
    return len(zlib.compress(s.encode("utf-8")))

def ncd(x, y):
    cx, cy, cxy = csize(x), csize(y), csize(x + y)
    return (cxy - min(cx, cy)) / max(cx, cy)

M = pd.DataFrame([[ncd(text[a], text[b]) for b in archs] for a in archs], index=archs, columns=archs)
M.to_csv(C.PROC / "ncd_matrix.csv")

pairs = []
for i in range(1, 9):
    for var in (f"NS{i}", f"SS{i}"):
        p, v = seqs[f"S{i}"], seqs[var]
        swapped = [k + 1 for k, (x, y) in enumerate(zip(p, v)) if x != y] if var.startswith("NS") else []
        pairs.append({"parent": f"S{i}", "variant": var, "type": var[:2], "NCD": M.loc[f"S{i}", var],
                      "swapped_positions": ",".join(map(str, swapped)) if swapped else "",
                      "variant_equals_parent_without_first_layer": var.startswith("SS") and v == p[1:]})
pv = pd.DataFrame(pairs)
pv.to_csv(C.PROC / "parent_variant_ncd.csv", index=False)

# ---------------------------------------------------------------- directed transitions
# Forward  = polycation deposited on a polyanion   (A->X, A->Y, B->X, B->Y)
# Reverse  = polyanion deposited on a polycation   (X->A, Y->A, X->B, Y->B)
FWD, REV = ["AX", "AY", "BX", "BY"], ["XA", "YA", "XB", "YB"]
TR = FWD + REV
LAB = {t: f"{C.NAME[t[0]]}->{C.NAME[t[1]]}" for t in TR}

counts = pd.DataFrame(
    [{t: sum(1 for k in range(len(seqs[a]) - 1) if seqs[a][k] + seqs[a][k + 1] == t) for t in TR} for a in archs],
    index=archs)
assert (counts.sum(axis=1) == comp.set_index("architecture").n_layers - 1).all(), "non-alternating sequence"
frac = counts.div(counts.sum(axis=1), axis=0)
counts.rename(columns=LAB).rename_axis("architecture").to_csv(C.PROC / "transition_counts.csv")
frac.rename(columns=LAB).rename_axis("architecture").to_csv(C.PROC / "transition_fractions.csv")

# Pairing-bias descriptors: +1 = only ALG/CHI and PEC/GEL pairings, -1 = only ALG/GEL and PEC/CHI
bias = pd.DataFrame({
    "B_A_to_C": ((counts.AX + counts.BY) - (counts.AY + counts.BX)) / counts[FWD].sum(axis=1),
    "B_C_to_A": ((counts.XA + counts.YB) - (counts.YA + counts.XB)) / counts[REV].sum(axis=1),
}).rename_axis("architecture")
bias.to_csv(C.PROC / "pairing_bias.csv")

# D_pair = total-variation distance between transition-fraction vectors
dp = []
for i in range(1, 9):
    for var in (f"NS{i}", f"SS{i}"):
        d = counts.loc[var] - counts.loc[f"S{i}"]
        dp.append({"parent": f"S{i}", "variant": var, "type": var[:2],
                   "D_pair": 0.5 * (frac.loc[f"S{i}"] - frac.loc[var]).abs().sum(),
                   "changed_counts": ", ".join(f"{LAB[t]} {int(v):+d}" for t, v in d.items() if v != 0)})
pd.DataFrame(dp).to_csv(C.PROC / "parent_variant_dpair.csv", index=False)

print(comp.groupby("family").n_layers.agg(["min", "max"]))
print(pv[["parent", "variant", "NCD", "swapped_positions"]].round(4).to_string(index=False))
print("Mean NCD  S->NS: %.4f   S->SS: %.4f" % (pv[pv.type == "NS"].NCD.mean(), pv[pv.type == "SS"].NCD.mean()))
