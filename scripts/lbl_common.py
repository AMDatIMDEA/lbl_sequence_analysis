"""
lbl_common.py
=============
Shared paths, settings and helper functions used by every analysis script.

Nothing in this file writes to disk; the numbered scripts (01_..., 02_..., ...)
import from here and produce the files in ``processed_data/`` and ``figures/``.

Conventions used throughout the repository
------------------------------------------
Polymer codes (as used in the deposition program and in ``raw_data/sequences``):
    A = alginate (ALG, polyanion)       B = pectin (PEC, polyanion)
    X = chitosan (CHI, polycation)      Y = gelatin (GEL, polycation)

Architecture families used in the paper:
    P   binary references (P1 Alg/Chi, P2 Pec/Gel, P3 Alg/Gel, P4 Pec/Chi)
    S   16 four-component sequences (S1-S16)
    NS  near-neighbour variants of S1-S8 (one polyanion-polyanion swap)
    SS  shifted variants of S1-S8 (first deposited layer removed)
Families that were fabricated but are NOT used in the paper (R/I, NR/NI) are
kept only in ``archive/``.

Tensile curves: strain is stored in % and stress in MPa (as exported by the DMA).
Therefore the quantities computed exactly as in the paper have these units:
    slope of the stress-strain fit  -> MPa per % strain  (multiply by 100 for MPa)
    area under the curve            -> MPa * %           (divide by 100 for MJ m^-3)
Both the paper-convention values and the SI-unit conversions are written to
``processed_data``.
"""

from pathlib import Path
import ast

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Paths (all relative to the repository root, so the repo can be moved)
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw_data"
PROC = ROOT / "processed_data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

# processed_data is split by topic (one sub-folder per type of result)
PROC_THICK = PROC / "thickness"          # 01
PROC_MECH = PROC / "mech_properties"     # 02
PROC_SEQ = PROC / "sequences"            # 03: composition, NCD, transitions
PROC_PV = PROC / "parent_variant"        # 03: parent-variant NCD and D_pair
PROC_PAIR = PROC / "pair_statistics"     # 03 (pairing bias) and 08
PROC_STATS = PROC / "statistics"         # 04, 05
PROC_GROWTH = PROC / "growth_test"       # 06
PROC_REP = PROC / "repeatability_test"   # 07
for _d in (PROC_THICK, PROC_MECH, PROC_SEQ, PROC_PV, PROC_PAIR, PROC_STATS, PROC_GROWTH, PROC_REP):
    _d.mkdir(parents=True, exist_ok=True)

FAMILIES = ["P", "S", "NS", "SS"]
NAME = {"A": "ALG", "B": "PEC", "X": "CHI", "Y": "GEL"}

# Tensile curves of all films (sheet "curves": one row per recorded point)
TENSILE_FILE = RAW / "tensile_curves" / "tensile_curves.xlsx"

# --------------------------------------------------------------------------
# Mechanical properties from one stress-strain curve
# (identical to the original paper_analysis/mech_analysis.py)
# --------------------------------------------------------------------------
def young_slope(strain, stress, stress_min, stress_max):
    """Least-squares slope of stress vs strain for points with
    stress_min <= stress <= stress_max. Returns MPa per % strain."""
    mask = (stress >= stress_min) & (stress <= stress_max)
    if mask.sum() < 5:
        return np.nan
    x, y = strain[mask], stress[mask]
    order = np.argsort(x)
    return np.polyfit(x[order], y[order], 1)[0]


def curve_properties(strain, stress, window):
    """All properties of one curve, in the paper convention."""
    return {
        "E_slope_MPa_per_pct": young_slope(strain, stress, *window),
        "UTS_MPa": np.nanmax(stress),
        "strain_at_break_pct": np.nanmax(strain),
        "toughness_MPa_pct": np.trapezoid(stress, strain),
    }


# --------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------
def read_tensile_curves(path=TENSILE_FILE):
    """Read the tensile workbook (sheet "curves").

    Columns: Family, Architecture, Film, Point, Strain (%), Stress (MPa).
    Returns a dict {(architecture, film_number): (strain_array, stress_array)},
    with the points in the order in which they were recorded.
    """
    df = pd.read_excel(path, sheet_name="curves").sort_values(["Architecture", "Film", "Point"])
    return {(a, int(f)): (g["Strain (%)"].to_numpy(float), g["Stress (MPa)"].to_numpy(float))
            for (a, f), g in df.groupby(["Architecture", "Film"], sort=False)}


def read_sequences(path=RAW / "sequences" / "sequences.txt"):
    """{architecture: list of polymer codes in deposition order}."""
    seqs = {}
    for line in Path(path).read_text().splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            seqs[k.strip()] = ast.literal_eval(v.strip())
    return seqs


def arch_sort_key(name):
    fam = name.rstrip("0123456789")
    return (FAMILIES.index(fam) if fam in FAMILIES else 99, int(name[len(fam):]))
