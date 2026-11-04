# LbL sequence library — data and analysis

Data and analysis code for the manuscript on sequence-controlled Layer-by-Layer (LbL)
multilayer films of alginate (ALG), pectin (PEC), chitosan (CHI) and gelatin (GEL),
assembled with the RoboMAT automated dip-coating platform.

Everything in `processed_data/` and `figures/` is generated from `raw_data/` by the
scripts in `scripts/`. Running

```bash
cd scripts
python run_all.py
```

rebuilds all of it (Python ≥ 3.10 with numpy, pandas, scipy, matplotlib, openpyxl).
The film-level table produced by step 02 is **identical** to the `mechprop.xlsx` used for
the manuscript (all 108 films, all properties; checked to machine precision).

> **Status:** organised from the first part of the original files (upload part 1).
> Items marked *pending part 2* are expected in the second upload.

---

## 1. What is in the paper

| Set | Architectures | Films | Description |
|---|---|---|---|
| P  | P1–P4   | 3 each | Binary references: P1 Alg/Chi, P2 Pec/Gel, P3 Alg/Gel, P4 Pec/Chi |
| S  | S1–S16  | 3 each | Four-component sequences, 100 layers, 25 of each polymer |
| NS | NS1–NS8 | 3 each | Near-neighbour variants of S1–S8: two polyanion layers swapped |
| SS | SS1–SS8 | 3 each | Shifted variants of S1–S8: first deposited layer removed (99 layers) |
| Growth test | 4 binary systems | 1 per point | Thickness after 5, 10, 15, 25, 35 and 50 bilayers |
| Repeatability test | Alg/Chi | 10 | Ten independently prepared films |

**Not used in the paper:** the R/I sets (16 sequences) and their near-neighbour variants
NR/NI (8). Their data are not included in this repository (only in the original, unorganised files). In the original NS/SS thickness file they are the
labels N9–N16, which were left out of `thickness.xlsx`.

## 2. Folder structure

```
raw_data/                     data as measured (never modified by the scripts)
  sequences/sequences.txt       deposition order of the 36 architectures (codes A, B, X, Y; read by the scripts)
  sequences/sequences_chem.txt  the same sequences written with polymer names (Alg, Pec, Chi, Gel)
  thickness/thickness.xlsx      all 108 films: readings, film means, architecture means (µm)
  tensile_curves/tensile_curves.xlsx        stress-strain curves of all 108 films (one row per point)
  growth_test/growth_test_thickness.xlsx
  repeatability_test/repeatability_tensile_curves.xlsx, repeatability_thickness.xlsx
scripts/                      annotated analysis, numbered in execution order
processed_data/               outputs of the scripts (CSV unless stated), one sub-folder per topic
  thickness/                    film thickness (01)
  mech_properties/              tensile properties of films and architectures (02)
  sequences/                    composition, NCD matrix, transition counts/fractions (03)
  parent_variant/               NCD and D_pair of every S -> NS/SS pair (03)
  pair_statistics/              pairing-bias descriptors (03) and SI Section 9 tables (08)
  statistics/                   group statistics and parent-variant changes (04, 05)
  growth_test/                  growth data and fits (06)
  repeatability_test/           repeatability results (07)
figures/                      figures produced by the scripts
```

### Where each raw file came from (original name → new name)

| New file | Original file |
|---|---|
| `raw_data/sequences/sequences.txt` | `sequences_final.txt` as corrected in Oct 2026 (NS2 = swap 5↔21, confirmed by Alex as deposited) — see Known issues 3 |
| `raw_data/thickness/thickness.xlsx` | `paper_analysis/OLD_DAT/Thickness S1-S16 and pure new.xlsx` (P, S) + `original_data/THICKNESS 8X8X8.xlsx` (NS, SS), merged into one workbook (one row per reading; NR left out; values unchanged) |
| `raw_data/tensile_curves/tensile_curves.xlsx` | merged from `paper_analysis/OLD_DAT/P_NEW.xlsx` (P), `paper_analysis/OLD_DAT/S_DMA_NEW.xlsx` (S), `paper_analysis/NS_DMA.xlsx` (NS) and `paper_analysis/SS_DMA.xlsx` (SS); one row per point, values unchanged |
| `raw_data/growth_test/growth_test_thickness.xlsx` | `original_data/thickness rep..xlsx`, reorganised (one row per film, labelled sheets, chart removed; values unchanged) |
| `raw_data/repeatability_test/repeatability_tensile_curves.xlsx` | `new_analysis/rep_DMA.xlsx`, reorganised (one row per point, README sheet; values unchanged) |
| `raw_data/repeatability_test/repeatability_thickness.xlsx` | `new_analysis/rep_thickness.xlsx`, reorganised (one row per reading, per-film summary, README sheet; values unchanged) |

These are the only versions that reproduce the manuscript values. Several files with the
same or similar names exist among the original files and give **different** numbers, e.g.
`final_analysis/S_DMA_NEW.xlsx`, `Thickness S1-S16 and pure 1.xlsx` (different S13–S16 and P).

Label conventions inside the raw workbooks:
- Original NS/SS thickness file: **N1–N8 = NS1–NS8**, N9–N16 = NR1–NR8 (not used), **S1–S8 = SS1–SS8**
  (kept in the column "Original label" of `thickness.xlsx`).
- Sequence codes: A = ALG, B = PEC (polyanions); X = CHI, Y = GEL (polycations).

## 3. Scripts

| Script | Does | Main outputs |
|---|---|---|
| `lbl_common.py` | Paths, settings, property definitions, readers | – |
| `01_thickness.py` | Reads all thickness readings | `thickness/thickness_readings.csv`, `thickness/thickness_films.csv` |
| `02_mechanical_properties.py` | Properties of every curve, joined with thickness | `mech_properties/mechanical_properties_films.csv`, `mech_properties/mechanical_properties_architectures.csv`, `mech_properties/mechprop_paper_format.xlsx`, `figures/stress_strain_*.png` |
| `03_sequence_analysis.py` | Composition, NCD (Eq. 1), directed transitions, pairing bias, D_pair | `sequences/` (`sequence_composition.csv`, `ncd_matrix.csv`, `transition_counts.csv`, `transition_fractions.csv`), `parent_variant/` (`parent_variant_ncd.csv`, `parent_variant_dpair.csv`), `pair_statistics/pairing_bias.csv` |
| `04_statistics.py` | Descriptive statistics, Welch ANOVA, Tukey HSD (S vs SS vs NS), Pearson matrix | `statistics/*.csv`, `figures/FigS3_*.png`, `figures/FigS4_*.png` |
| `05_parent_variant.py` | S→NS and S→SS changes, Wilcoxon, Spearman vs NCD | `statistics/parent_variant_changes.csv`, `statistics/SS_wilcoxon.csv`, `statistics/NS_spearman_vs_NCD.csv` |
| `06_growth_test.py` | Growth curves and model fits | `growth_test/*.csv`, `figures/LbL_growth_models.png` |
| `07_repeatability_test.py` | Ten Alg/Chi films, CVs with/without film 7 | `repeatability_test/*.csv`, `figures/repeatability_curves.png` |
| `08_pair_statistics.py` / `.ipynb` | SI Section 9 (transitions heat map, pairing-bias regression, D_pair) | `pair_statistics/*.csv`, `figures/pair_stats_*.png/.pdf` |
| `run_all.py` | Runs 01–08 | everything |

### Property definitions (per stress–strain curve; strain in %, stress in MPa)

| Column | Definition | Paper column name |
|---|---|---|
| `thickness_um` | film mean recorded in the thickness workbook (= mean of its readings, except NS8 film 3) (µm) | `Thickness_mm` (values are in µm) |
| `E_slope_MPa_per_pct` | slope of a linear fit to the initial part of the curve (stress range set in `E_FIT_RANGE` at the top of `02_mechanical_properties.py`) | `Young_Modulus_MPa` |
| `E_MPa` | `E_slope × 100` — Young's modulus in MPa | – |
| `UTS_MPa` | maximum stress | `UTS_MPa` |
| `strain_at_break_pct` | maximum strain | `Elongation_%` |
| `toughness_MPa_pct` | trapezoidal area under the curve | `Toughness` |
| `toughness_MJ_m3` | `toughness_MPa_pct / 100` | – |

`mechprop_paper_format.xlsx` keeps the original column names so that older notebooks
(and `08_pair_statistics`) can read it unchanged.
