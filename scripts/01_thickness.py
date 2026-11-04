"""
01_thickness.py
===============
Reads the film-thickness workbook and writes tidy tables for the other scripts.

Input   raw_data/thickness/thickness.xlsx
            sheet "readings": one row per reading
                Family, Architecture, Film, Reading, Thickness (µm)
                (P and S films: 5 readings; NS and SS films: 3 readings)
            sheet "films": one row per film, "Thickness (µm)" = film thickness used in the manuscript
Outputs (processed_data/thickness/)
    thickness_readings.csv   one row per reading
    thickness_films.csv      one row per film: thickness_um (used downstream) and n_readings
The film thickness is the mean of its readings, except NS8 film 3, for which the manuscript
used 33.0 µm (mean of readings 32.67 µm). The script prints any film where the two differ.
"""
import pandas as pd
from lbl_common import RAW, PROC_THICK, arch_sort_key

FILE = RAW / "thickness" / "thickness.xlsx"

r = pd.read_excel(FILE, sheet_name="readings")
readings = pd.DataFrame({"architecture": r["Architecture"], "film": r["Film"], "reading": r["Reading"],
                         "thickness_um": r["Thickness (µm)"], "family": r["Family"]})
readings = readings.sort_values(["architecture", "film", "reading"],
                                key=lambda s: s.map(arch_sort_key) if s.name == "architecture" else s)

recorded = pd.read_excel(FILE, sheet_name="films")[["Architecture", "Film", "Thickness (µm)"]]
recorded.columns = ["architecture", "film", "thickness_um"]

films = (readings.groupby(["family", "architecture", "film"], sort=False)
         .agg(thickness_um_recomputed=("thickness_um", "mean"), n_readings=("thickness_um", "size"))
         .reset_index()
         .merge(recorded, on=["architecture", "film"], how="left", validate="one_to_one"))

diff = films[(films.thickness_um - films.thickness_um_recomputed).abs() > 1e-6]
if len(diff):
    print("Films whose thickness differs from the mean of their readings:")
    print(diff[["architecture", "film", "thickness_um", "thickness_um_recomputed"]].round(3).to_string(index=False))

readings.to_csv(PROC_THICK / "thickness_readings.csv", index=False)
films[["family", "architecture", "film", "thickness_um", "n_readings"]].to_csv(PROC_THICK / "thickness_films.csv", index=False)

print(films.groupby("family").agg(films=("film", "size"), readings_per_film=("n_readings", "max")))
print("Written thickness_readings.csv and thickness_films.csv")
