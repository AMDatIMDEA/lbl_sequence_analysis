"""
08_pair_statistics.py
=====================
Plain-script version of 08_pair_statistics.ipynb (identical code; use whichever you prefer).
Analysis of sequence dependence, Supporting Information Section 9:
  9.1 directed nearest-neighbour deposition transitions (heat map, coupling)
  9.2 pairing-bias descriptors B_A->C, B_C->A and regression on the S library
  9.3 parent-variant distance D_pair vs changes in mechanical properties
Run from the scripts/ folder after 01-03. Tables go to processed_data/pair_statistics/, figures to figures/.
Mechanical properties are read in the paper convention (mech_properties/mechprop_paper_format.xlsx).
"""
import matplotlib
matplotlib.use("Agg")
def display(x):
    print(x.to_string() if hasattr(x,'to_string') else x)
import ast, re, os
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy import stats

# ---------------- SETTINGS ----------------
SEQ_FILE  = "../raw_data/sequences/sequences.txt"
MECH_FILE = "../processed_data/mech_properties/mechprop_paper_format.xlsx"
OUT_DIR   = "../processed_data/pair_statistics"   # tables
FIG_DIR   = "../figures"                          # figures
AGG       = "mean"      # how the 3 replicates are combined: "mean" or "median"
DPI       = 600         # resolution of the PNG figures

# column in mechprop.xlsx -> label used in tables and figures
PROPS = {"Young_Modulus_MPa": "E (MPa)",
         "UTS_MPa":           "UTS (MPa)",
         "Elongation_%":      "ε (%)",
         "Toughness":         "Toughness"}
THICK = "Thickness_mm"

# polymer codes used in the sequence file
NAME    = {"A": "ALG", "B": "PEC", "X": "CHI", "Y": "GEL"}
ANIONS  = ("A", "B")
CATIONS = ("X", "Y")
# the 8 directed transitions
FWD = ["AX", "AY", "BX", "BY"]          # polyanion -> polycation
REV = ["XA", "YA", "XB", "YB"]          # polycation -> polyanion
TRANS = FWD + REV
LABEL = {t: f"{NAME[t[0]]}→{NAME[t[1]]}" for t in TRANS}

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
display(pd.set_option("display.width", 200, "display.max_columns", 40))

# ---------------- FIGURE STYLE (journal-ready) ----------------
FONT_PANEL     = 20     # base font size (pt) of the individual panels
FONT_COMPOSITE = 18      # base font size (pt) of the composite figures (≈ 183 mm wide)

def set_style(base):
    # all font sizes are derived from one base size, so a single number controls the whole figure
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "xtick.labelsize": base - 1, "ytick.labelsize": base - 1, "legend.fontsize": base - 1,
        "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,      # keep text editable in Illustrator/Inkscape
        "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
    })

def fs(delta=0):
    # font size relative to the current base size
    return mpl.rcParams["font.size"] + delta

set_style(FONT_PANEL)
SEQ_CMAP  = "Blues"                 # single-hue sequential ramp (magnitude)
C_NS, C_SS = "#2a78d6", "#eb6834"   # categorical: NS = blue, SS = orange (CVD-validated pair)
INK, MUTED = "#1a1a1a", "#6b6b6b"

def label_points(ax, xs, ys, texts, side="right", marker_size=40):
    # Direct labels without overlaps: points sharing (almost) the same x are labelled in a column
    # on one side, pushed apart vertically so that labels never touch. Call after setting limits.
    fontsize = fs(-2)
    bb = ax.get_position(); W, H = ax.figure.get_size_inches()
    ax_w, ax_h = bb.width * W * 72, bb.height * H * 72                 # axes size in points
    gap = 1.15 * fontsize / ax_h                                       # one text line, in axes fraction
    dx = (np.sqrt(marker_size) / 2 + 2) / ax_w                         # clear the marker radius
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    fx = [(x - x0) / (x1 - x0) for x in xs]; fy = [(y - y0) / (y1 - y0) for y in ys]
    groups = {}
    for i, x in enumerate(fx):
        groups.setdefault(round(x / 0.02), []).append(i)
    for idx in groups.values():
        idx = sorted(idx, key=lambda i: fy[i]); pos = []
        for i in idx:
            pos.append(max(fy[i], pos[-1] + gap) if pos else fy[i])
        overflow = max(0, pos[-1] - 0.98)
        pos = [q - overflow for q in pos]
        # side="auto": label the left-most column on its left, every other column on its right
        sd = side if side != "auto" else ("left" if min(fx[i] for i in idx) < 0.2 else "right")
        for i, q in zip(idx, pos):
            xt = fx[i] + (dx if sd == "right" else -dx)
            ax.text(xt, q, texts[i], transform=ax.transAxes,
                    ha="left" if sd == "right" else "right", va="center", fontsize=fontsize, color=MUTED)
            if abs(q - fy[i]) > 0.5 * gap:          # label was pushed away: thin leader line to its point
                x_edge = fx[i] + (0.6 if sd == "right" else -0.6) * dx
                ax.plot([x_edge, xt + (-0.15 if sd == "right" else 0.15) * dx], [fy[i], q], transform=ax.transAxes,
                        color=MUTED, lw=0.4, alpha=0.7, clip_on=False, zorder=1)

def save(fig, name):
    # save every figure as vector PDF + high-resolution PNG
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=DPI)

def panel_heatmap(ax, ids, orientation="vertical", cbar="bottom"):
    # (a) colour = transition fraction p_i, number = raw count
    # orientation="vertical": sequences as rows; "horizontal": sequences as columns (wide panel)
    M = P.loc[ids, TRANS].values.astype(float)
    C = N.loc[ids, TRANS].values
    if orientation == "horizontal":
        M, C = M.T, C.T
    im = ax.imshow(M, cmap=SEQ_CMAP, aspect="auto", vmin=M.min() * 0.6, vmax=M.max())
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, int(C[i, j]), ha="center", va="center", fontsize=fs(-2.5),
                    color="white" if M[i, j] > 0.6 * M.max() else INK)
    tlabels, slabels = [LABEL[t] for t in TRANS], list(ids)
    ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    if orientation == "vertical":
        ax.set_xticks(range(8), tlabels, rotation=60, ha="right", rotation_mode="anchor")
        ax.set_yticks(range(len(ids)), slabels)
        ax.axvline(3.5, color="white", lw=2.5)                        # separate the two directions
        for x, txt in ((1.5, "polyanion →\npolycation"), (5.5, "polycation →\npolyanion")):
            ax.text(x, -0.7, txt, ha="center", va="bottom", fontsize=fs(-1.5), color=MUTED)
    else:
        ax.set_yticks(range(8), tlabels)
        ax.set_xticks(range(len(ids)), slabels, rotation=90)
        ax.axhline(3.5, color="white", lw=2.5)
        # group labels on the left, horizontal, in two lines (they always fit in the height of 4 rows)
        for y, txt in ((1.5, "polyanion →\npolycation"), (5.5, "polycation →\npolyanion")):
            ax.annotate(txt, xy=(0, y), xycoords=("axes fraction", "data"), xytext=(-fs(1) * 6.2, 0),
                        textcoords="offset points", ha="right", va="center", fontsize=fs(-5), color=MUTED)
    if cbar == "bottom":
        cb = plt.colorbar(im, ax=ax, orientation="horizontal", fraction=0.04, pad=0.12 if orientation == "vertical" else 0.25, aspect=35)
    else:
        cb = plt.colorbar(im, ax=ax, orientation="vertical", fraction=0.03, pad=0.02, aspect=20)
    cb.set_label("transition fraction $p_i$"); cb.outline.set_linewidth(0.5); cb.ax.tick_params(labelsize=fs(-2))
    return im

def _is_bottom(n, grid): return n // grid[1] == grid[0] - 1
def _is_left(n, grid):   return n % grid[1] == 0

def panel_bias(axes, grid=(2, 2), msize=55, labels=True):
    # (b) S1-S16 in (B_fwd, B_rev) space, one subplot per property, coloured by its value
    xb, yb = B.loc[S_ids, "B_fwd"], B.loc[S_ids, "B_rev"]
    pad = 0.06
    for n, ((prop, lab), ax) in enumerate(zip(PROPS.items(), axes)):
        v = Y.loc[S_ids, prop]
        sc = ax.scatter(xb, yb, c=v, cmap=SEQ_CMAP, s=msize, edgecolor=INK, linewidth=0.5,
                        vmin=v.min() - 0.15 * (v.max() - v.min()), vmax=v.max(), zorder=2)
        r = reg.loc[lab]
        ax.set_title(lab, fontsize=fs(-0.5), pad=fs(1) * 3.4)   # title well above the subtitle
        ax.text(0.5, 1.04, f"$Q^2_{{LOO}}$ = {r['Q2 LOO']:.2f}   adj. $R^2$ = {r['adj R2']:.2f}", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=fs(-2.5), color=MUTED)
        ax.set_xlim(xb.min() - 2.5 * pad, xb.max() + 2 * pad); ax.set_ylim(yb.min() - pad, yb.max() + pad)
        ax.axhline(0, color="#d0d0d0", lw=0.6, zorder=0); ax.axvline(0, color="#d0d0d0", lw=0.6, zorder=0)
        ax.set_xlabel("$B_{A\\to C}$ (forward)" if _is_bottom(n, grid) else "")
        ax.set_ylabel("$B_{C\\to A}$ (reverse)" if _is_left(n, grid) else "")
        if not _is_bottom(n, grid): ax.tick_params(labelbottom=False)
        cb = plt.colorbar(sc, ax=ax, fraction=0.06, pad=0.03)
        cb.outline.set_linewidth(0.5); cb.ax.tick_params(labelsize=fs(-2))
        # one label per distinct point; sequences with identical (B_fwd, B_rev) share a label, e.g. "8,10"
        if not labels: continue
        coords = {}
        for k in S_ids:
            coords.setdefault((round(xb[k], 4), round(yb[k], 4)), []).append(str(num(k)))
        label_points(ax, [c[0] for c in coords], [c[1] for c in coords],
                     [",".join(v) for v in coords.values()], side="auto", marker_size=msize)

def panel_dpair(axes, grid=(2, 2), msize=45, legend=True, label_min=0.0):
    # label_min: only points above this fraction of the y-range get a sequence number (0 = label all)
    # (c) pair-statistics distance vs |change| in each property, for S->NS and S->SS
    xmax = pairs.D_pair.max() * 1.2
    for n, ((prop, lab), ax) in enumerate(zip(PROPS.items(), axes)):
        for fam, color, marker in (("NS", C_NS, "o"), ("SS", C_SS, "o")):
            d = pairs[pairs.type == fam]
            ax.errorbar(d.D_pair, d[f"|Δ| {prop}"], yerr=d[f"SE {prop}"], fmt="none",
                        ecolor=color, elinewidth=0.6, capsize=1.5, alpha=0.6, zorder=1)
            ax.scatter(d.D_pair, d[f"|Δ| {prop}"], color=color, marker=marker, s=msize,
                       edgecolor="white", linewidth=0.7, zorder=2, label=f"S → {fam}")
        ax.set_ylim(0, (pairs[f"|Δ| {prop}"] + pairs[f"SE {prop}"]).max() * 1.05)
        ax.set_xlim(-0.2 * xmax, xmax)
        ax.set_ylabel(f"|Δ {lab}|")
        ax.set_xlabel("$D_{pair}$" if _is_bottom(n, grid) else "")
        if not _is_bottom(n, grid): ax.tick_params(labelbottom=False)
        # number = parent sequence; SS labels to the left of their column, NS labels to the right
        # number = parent sequence; SS labels to the right of their column; NS labels to the right,
        # except the left-most NS column (e.g. D_pair = 0), labelled on its left
        ymax = ax.get_ylim()[1]
        for fam, side in (("SS", "right"), ("NS", "auto")):
            d = pairs[(pairs.type == fam) & (pairs[f"|Δ| {prop}"] >= label_min * ymax)]
            label_points(ax, d.D_pair.tolist(), d[f"|Δ| {prop}"].tolist(),
                         [r.lstrip("S") for r in d.parent], side=side, marker_size=msize)
    if legend:
        # legend above the first subplot, outside the data area, so it never hides a point
        axes[0].legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2,
                       handletextpad=0.2, columnspacing=1.0, borderaxespad=0)

def add_letters(fig, items, dx=-0.01, dy=0.005, align=True):
    # bold panel letters at the top of each panel's bounding box (incl. titles and labels);
    # align=True puts all letters on the same vertical line (the left-most edge of all panels)
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    boxes = {letter: [ax.get_tightbbox(r).transformed(inv) for ax in axes] for letter, axes in items}
    x_common = min(b.x0 for bs in boxes.values() for b in bs)
    for letter, bs in boxes.items():
        x0 = x_common if align else min(b.x0 for b in bs)
        y1 = max(b.y1 for b in bs)
        fig.text(x0 + dx, y1 + dy, letter, fontsize=fs(5), fontweight="bold", ha="left", va="bottom")

def load_sequences(path):
    # each line of the file looks like:  S1: ['B', 'X', 'A', ...]
    seqs = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*(\w+)\s*:\s*(\[.*\])", line)
            if m:
                seqs[m.group(1)] = list(ast.literal_eval(m.group(2)))
    return seqs

seqs   = load_sequences(SEQ_FILE)
family = lambda k: re.match(r"[A-Z]+", k).group()             # "S12" -> "S"
num    = lambda k: int(re.sub(r"\D", "", k))                  # "S12" -> 12
order  = lambda k: ({"P": 0, "S": 1, "NS": 2, "SS": 3}[family(k)], num(k))
S_ids  = sorted([k for k in seqs if family(k) == "S"], key=order)
NS_ids = sorted([k for k in seqs if family(k) == "NS"], key=order)
SS_ids = sorted([k for k in seqs if family(k) == "SS"], key=order)

# mechanical data: one row per replicate -> one value per sequence
mech = pd.read_excel(MECH_FILE)
cols = list(PROPS) + [THICK]
Y   = mech.groupby("Sample")[cols].agg(AGG)     # mean (or median) of the 3 replicates
Ysd = mech.groupby("Sample")[cols].std()        # standard deviation of the replicates
Yn  = mech.groupby("Sample")[cols].count()      # number of replicates

# number of layers and of each polymer, per family
comp = pd.DataFrame({k: {"layers": len(s), **{NAME[p]: s.count(p) for p in "ABXY"}} for k, s in seqs.items()}).T
display(comp.groupby(comp.index.map(family)).agg(["min", "max"]).loc[["S", "NS", "SS"]])

def transition_counts(seq):
    c = Counter(seq[i] + seq[i + 1] for i in range(len(seq) - 1))
    return pd.Series({t: c.get(t, 0) for t in TRANS})

N = pd.DataFrame({k: transition_counts(s) for k, s in seqs.items()}).T.loc[sorted(seqs, key=order)]   # counts
P = N.div(N.sum(axis=1), axis=0)                                                                      # fractions

table = N.rename(columns=LABEL)
table.insert(0, "Layers", [len(seqs[k]) for k in table.index])
table["Total"] = N.sum(axis=1)
table.to_csv(os.path.join(OUT_DIR, "transition_counts.csv"))

# LaTeX version for the SI (S, NS, SS)
def latex_counts():
    L = [r"\begin{table}[ht]", r"\centering", r"\small",
         r"\caption{Directed nearest-neighbour deposition transitions between consecutive layers (deposition order, substrate first). S and NS: 100 layers (99 transitions); SS: 99 layers (98 transitions).}",
         r"\label{tab:transitions}", r"\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{lcccccccc}", r"\toprule",
         r" & \multicolumn{4}{c}{Polyanion $\rightarrow$ polycation} & \multicolumn{4}{c}{Polycation $\rightarrow$ polyanion} \\",
         r"\cmidrule(lr){2-5} \cmidrule(lr){6-9}",
         "Sequence & " + " & ".join(LABEL[t].replace("→", r"$\rightarrow$") for t in TRANS) + r" \\", r"\midrule"]
    prev = None
    for k in S_ids + NS_ids + SS_ids:
        if prev and family(k) != prev: L.append(r"\midrule")
        prev = family(k)
        L.append(f"{k} & " + " & ".join(str(int(N.loc[k, t])) for t in TRANS) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(L)
open(os.path.join(OUT_DIR, "transition_counts.tex"), "w").write(latex_counts())

display(table)

Ns = N.loc[S_ids].astype(float)
same = [(LABEL[a], LABEL[b]) for i, a in enumerate(TRANS) for b in TRANS[i + 1:] if (Ns[a] == Ns[b]).all()]
print("Columns that are identical in S1–S16:", same)
print("ALG→CHI + ALG→GEL in S1–S16:", sorted(int(v) for v in (Ns.AX + Ns.AY).unique()))
display(Ns.rename(columns=LABEL).corr().round(2))

B = pd.DataFrame({
    "B_fwd": ((N.AX + N.BY) - (N.AY + N.BX)) / N[FWD].sum(axis=1),
    "B_rev": ((N.XA + N.YB) - (N.YA + N.XB)) / N[REV].sum(axis=1),
})

# how much of each count do the two descriptors reproduce? (linear fit, S1-S16)
X1 = np.column_stack([np.ones(len(S_ids)), B.loc[S_ids].values])
r2 = {}
for t in TRANS:
    y = Ns[t].values
    res = y - X1 @ np.linalg.lstsq(X1, y, rcond=None)[0]
    r2[LABEL[t]] = 1 - (res ** 2).sum() / ((y - y.mean()) ** 2).sum()
print("R² of each count reproduced by (B_fwd, B_rev):")
print(pd.Series(r2).round(3).to_string())
display(B.loc[S_ids].round(3).T)

set_style(FONT_PANEL)
fig, ax = plt.subplots(figsize=(15.0, 4.0))
panel_heatmap(ax, S_ids, orientation="horizontal", cbar="right")
save(fig, "pair_stats_a_heatmap"); plt.show()

def ols(X, y):
    # linear regression with intercept; returns R2, adjusted R2, F-test p-value, leave-one-out Q2 and predictions
    X1 = np.column_stack([np.ones(len(y)), np.asarray(X, float).reshape(len(y), -1)])
    n, p = X1.shape
    beta = np.linalg.lstsq(X1, y, rcond=None)[0]
    res  = y - X1 @ beta
    h    = np.diag(X1 @ np.linalg.pinv(X1))            # leverage of each point
    loo  = res / (1 - h)                               # exact leave-one-out residuals
    sst  = ((y - y.mean()) ** 2).sum()
    r2   = 1 - (res ** 2).sum() / sst
    F    = (r2 / (p - 1)) / ((1 - r2) / (n - p))
    return dict(R2=r2, adjR2=1 - (1 - r2) * (n - 1) / (n - p), p_F=stats.f.sf(F, p - 1, n - p),
                Q2=1 - (loo ** 2).sum() / sst, pred_loo=y - loo, beta=beta)

D = B.loc[S_ids].join(Y.loc[S_ids])
rows = []
for prop, lab in PROPS.items():
    y  = D[prop].values
    m  = ols(D[["B_fwd", "B_rev"]].values, y)
    mt = ols(D[["B_fwd", "B_rev", THICK]].values, y)
    rf, pf = stats.spearmanr(D.B_fwd, y)
    rr, pr = stats.spearmanr(D.B_rev, y)
    rows.append({"Property": lab, "R2": m["R2"], "adj R2": m["adjR2"], "p (F)": m["p_F"], "Q2 LOO": m["Q2"],
                 "adj R2 (+thickness)": mt["adjR2"], "Q2 LOO (+thickness)": mt["Q2"],
                 "rho B_fwd": rf, "p rho B_fwd": pf, "rho B_rev": rr, "p rho B_rev": pr})
reg = pd.DataFrame(rows).set_index("Property")
reg.to_csv(os.path.join(OUT_DIR, "regression.csv"))
display(reg.round(3))

best_q2, best_adj = reg["Q2 LOO"].max(), reg["adj R2"].max()
any_signif = (reg["p (F)"] < 0.05).any()
if best_q2 > 0.5:
    verdict = "MOST: the pair descriptors explain most of the property variation."
elif any_signif or best_q2 > 0:
    verdict = "SOME: pair statistics contribute, but do not determine the mechanical state."
else:
    verdict = "LITTLE: the pair descriptors explain little; behaviour is not reducible to nearest-neighbour statistics."
print(f"adjusted R² from {reg['adj R2'].min():.2f} to {best_adj:.2f};  Q² LOO from {reg['Q2 LOO'].min():.2f} to {best_q2:.2f};  "
      f"smallest p (F) = {reg['p (F)'].min():.2f}")
print("Verdict ->", verdict)

set_style(FONT_PANEL)
fig, axes = plt.subplots(1, 4, figsize=(25.0, 6.0))
panel_bias(axes.ravel(), grid=(1, 4), msize=150)
fig.tight_layout(w_pad=2.0); save(fig, "pair_stats_b_pair_bias"); plt.show()

rows = []
for k in NS_ids + SS_ids:
    parent = f"S{num(k)}"
    rec = {"parent": parent, "variant": k, "type": family(k),
           "D_pair": 0.5 * (P.loc[parent] - P.loc[k]).abs().sum(),
           "changed counts": ", ".join(f"{LABEL[t]} {int(d):+d}" for t, d in (N.loc[k] - N.loc[parent]).items() if d)}
    for prop in PROPS:
        rec[f"|Δ| {prop}"] = abs(Y.loc[k, prop] - Y.loc[parent, prop])
        rec[f"SE {prop}"]  = np.sqrt(Ysd.loc[k, prop] ** 2 / Yn.loc[k, prop] + Ysd.loc[parent, prop] ** 2 / Yn.loc[parent, prop])
    rows.append(rec)
pairs = pd.DataFrame(rows)
pairs.to_csv(os.path.join(OUT_DIR, "parent_variant.csv"), index=False)

# does a larger change in pair statistics give a larger change in the property?
for prop, lab in PROPS.items():
    rho, p = stats.spearmanr(pairs.D_pair, pairs[f"|Δ| {prop}"])
    print(f"{lab:10s}  Spearman(D_pair, |Δ|) = {rho:+.2f}  (p = {p:.2f})")
display(pairs[["parent", "variant", "D_pair", "changed counts"] + [f"|Δ| {p}" for p in PROPS]].round(3))

set_style(FONT_PANEL)
fig, axes = plt.subplots(1, 4, figsize=(25.0, 5.0))
panel_dpair(axes.ravel(), grid=(1, 4), msize=130)
fig.tight_layout(w_pad=2.0); save(fig, "pair_stats_c_Dpair"); plt.show()

ss = pairs[pairs.type == "SS"]
summary = pd.DataFrame({lab: {"smallest |Δ|": ss[f"|Δ| {prop}"].min(), "in": ss.loc[ss[f"|Δ| {prop}"].idxmin(), "variant"],
                              "largest |Δ|": ss[f"|Δ| {prop}"].max(),  "in ": ss.loc[ss[f"|Δ| {prop}"].idxmax(), "variant"]}
                        for prop, lab in PROPS.items()}).T
print(f"D_pair of the SS variants: {ss.D_pair.min():.4f} – {ss.D_pair.max():.4f}")
display(summary)

rng = np.random.default_rng(42)

def perm_p(X, y, n_perm=10000):
    X1 = np.column_stack([np.ones(len(y)), X]); H = X1 @ np.linalg.pinv(X1)
    r2 = lambda v: 1 - ((v - H @ v) ** 2).sum() / ((v - v.mean()) ** 2).sum()
    null = np.array([r2(rng.permutation(y)) for _ in range(n_perm)])
    return (1 + (null >= r2(y)).sum()) / (n_perm + 1)

def holm(p):
    p = np.asarray(p); idx = np.argsort(p); out = np.empty(len(p)); run = 0
    for rank, i in enumerate(idx):
        run = max(run, (len(p) - rank) * p[i]); out[i] = min(1, run)
    return out

reg["p (perm)"] = [perm_p(D[["B_fwd", "B_rev"]].values, D[p].values) for p in PROPS]
reg["p (perm, Holm)"] = holm(reg["p (perm)"])
display(reg[["R2", "adj R2", "p (perm)", "p (perm, Holm)"]].round(3))

set_style(FONT_PANEL)
fig, axes = plt.subplots(1, 4, figsize=(25.0, 5.5))
for ax, (prop, lab) in zip(axes, PROPS.items()):
    y = D[prop].values; pred = ols(D[["B_fwd", "B_rev"]].values, y)["pred_loo"]
    ax.scatter(y, pred, s=130, color=C_NS, edgecolor="white", zorder=2)
    lo, hi = min(y.min(), pred.min()), max(y.max(), pred.max())
    ax.plot([lo, hi], [lo, hi], color="#999999", lw=1, ls="--", zorder=1)
    ax.axhline(y.mean(), color="#d0d0d0", lw=0.8, zorder=0)
    ax.set(title=f"{lab}   $Q^2$ = {reg.loc[lab, 'Q2 LOO']:.2f}", xlabel="measured", ylabel="predicted (LOO)")
fig.tight_layout(w_pad=2.0); save(fig, "pair_stats_extra_loo_predictions"); plt.show()

Y_median = mech.groupby("Sample")[cols].median()
rows = []
for name, Ymat, ids in (("S, mean", Y, S_ids), ("S, median", Y_median, S_ids), ("S+NS+SS, mean", Y, S_ids + NS_ids + SS_ids)):
    d = B.loc[ids].join(Ymat.loc[ids])
    for prop, lab in PROPS.items():
        m = ols(d[["B_fwd", "B_rev"]].values, d[prop].values)
        rows.append({"analysis": name, "n": len(ids), "property": lab, "adj R2": m["adjR2"], "Q2 LOO": m["Q2"], "p (F)": m["p_F"]})
robust = pd.DataFrame(rows)
display(robust.pivot(index="analysis", columns="property", values="adj R2").round(2))

def switching(seq):
    an = [c for c in seq if c in ANIONS]; ca = [c for c in seq if c in CATIONS]
    return {"switch_an": np.mean([an[i] != an[i + 1] for i in range(len(an) - 1)]),
            "switch_ca": np.mean([ca[i] != ca[i + 1] for i in range(len(ca) - 1)])}

SW = pd.DataFrame({k: switching(seqs[k]) for k in S_ids}).T
D2 = SW.join(Y.loc[S_ids])
rows = []
for prop, lab in PROPS.items():
    y = D2[prop].values; m = ols(D2[["switch_an", "switch_ca"]].values, y)
    rows.append({"Property": lab, "adj R2": m["adjR2"], "Q2 LOO": m["Q2"],
                 "p (perm)": perm_p(D2[["switch_an", "switch_ca"]].values, y),
                 "rho switch_an": stats.spearmanr(D2.switch_an, y)[0], "rho switch_ca": stats.spearmanr(D2.switch_ca, y)[0]})
sw = pd.DataFrame(rows).set_index("Property")
sw["p (perm, Holm)"] = holm(sw["p (perm)"])
sw["adj R2 (pair descriptors)"] = reg["adj R2"]          # side-by-side comparison with Part A
display(sw.round(3))

for d in (OUT_DIR, FIG_DIR):
    print("Files saved in", os.path.abspath(d))
    for f in sorted(os.listdir(d)): print("  ", f)

