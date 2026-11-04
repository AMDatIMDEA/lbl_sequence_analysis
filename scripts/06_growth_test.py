"""
06_growth_test.py
=================
Film growth of the four binary (pure) systems: thickness vs number of deposited layers.

Input   raw_data/growth_test/growth_test_thickness.xlsx, sheet "growth_data"
        One row per film: System, Code, Bilayers, Layers, Replicate, Thickness (mm), Thickness (µm).
        AX = Alg/Chi, AY = Alg/Gel, BX = Pec/Chi, BY = Pec/Gel.
        To add replicates, add rows with Replicate = 2, 3, ... (nothing to change here):
        the points become the mean of the replicates and the error bars their SD.
Outputs processed_data/growth_test/growth_data.csv   mean, SD and n per system and point (incl. origin)
        processed_data/growth_test/growth_fits.csv   fitted model, parameters and R^2 per system
        figures/LbL_growth_models.png                2x2 figure (600 dpi)

Growth models, fitted to the mean values with x = layers / 100:
    linear      T = m x + b               (Alg/Chi, Alg/Gel)
    saturating  T = a (1 - exp(-b x))     (Pec/Gel, Pec/Chi)
1 bilayer = 2 layers, so 100 layers = 50 bilayers. Set X_AXIS = "bilayers" to plot 0-50 bilayers.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import lbl_common as C

OUT = C.PROC / "growth_test"
OUT.mkdir(exist_ok=True)

# ==========================================================
# Settings
# ==========================================================

EXCEL_FILE = C.RAW / "growth_test" / "growth_test_thickness.xlsx"
FIGURE = C.FIG / "LbL_growth_models.png"

X_AXIS = "layers"            # "layers" (0-100) or "bilayers" (0-50)

# panel order: code in the Excel -> panel title
systems = {
    "AX": "Alg/Chi",
    "BY": "Pec/Gel",
    "AY": "Alg/Gel",
    "BX": "Pec/Chi",
}
colors = ["cornflowerblue", "darkorange", "seagreen", "firebrick"]

# ==========================================================
# Read the data
# ==========================================================

df = pd.read_excel(EXCEL_FILE, sheet_name="growth_data")
df["Thickness_um"] = df["Thickness (mm)"] * 1000      # recalculated here, not read from the formula column
df["x"] = df["Layers"] if X_AXIS == "layers" else df["Bilayers"]
X_MAX = 100 if X_AXIS == "layers" else 50

# mean, standard deviation and number of replicates per system and point
stats = (df.groupby(["Code", "x"])["Thickness_um"]
           .agg(["mean", "std", "count"])
           .reset_index())
stats["std"] = stats["std"].fillna(0)                  # one replicate -> no error bar

# ==========================================================
# Models
# ==========================================================

def linear(n, m, b):
    return m * n + b

def exp_safe(n, a, b):
    return a * (np.exp(b * n) - 1)

def sat_safe(n, a, b):
    return a * (1 - np.exp(-b * n))

model_map = {
    "AX": (linear, [35, 0]),
    "BY": (sat_safe, [30, 3]),
    "AY": (linear, [0.4, 0]),
    "BX": (sat_safe, [35, 2]),
}

# ==========================================================
# Style
# ==========================================================

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
})

fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
axs = axs.flatten()

# ==========================================================
# Fit and plot
# ==========================================================

table, fits = [], []

for i, (code, name) in enumerate(systems.items()):

    d = stats[stats["Code"] == code].sort_values("x")

    # add the origin (0 layers -> 0 µm)
    x = np.concatenate([[0], d["x"].values])
    y = np.concatenate([[0], d["mean"].values])
    sd = np.concatenate([[0], d["std"].values])
    n_rep = np.concatenate([[0], d["count"].values])

    # fit on the mean values (x normalised by its maximum)
    model, p0 = model_map[code]
    popt, _ = curve_fit(model, x / X_MAX, y, p0=p0, maxfev=20000)
    pred = model(x / X_MAX, *popt)
    r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)

    x_plot = np.linspace(0, X_MAX, 400)
    y_fit = model(x_plot / X_MAX, *popt)

    axs[i].errorbar(x, y, yerr=sd, fmt="o", ms=5, color=colors[i],
                    mec="black", mew=0.3, ecolor=colors[i], elinewidth=1,
                    capsize=3, zorder=3)
    axs[i].plot(x_plot, y_fit, "--", color=colors[i], linewidth=1.6)

    axs[i].set_title(name, pad=6)
    axs[i].set_xlim(0, X_MAX * 1.1)
    axs[i].set_ylim(0, 50)
    axs[i].set_xticks(np.linspace(0, X_MAX, 6))
    axs[i].grid(alpha=0.15)

    # panel letter (not bold)
    axs[i].text(0.03, 0.96, chr(65 + i), transform=axs[i].transAxes,
                fontsize=12, va="top")

    table += [{"system": name, "code": code, X_AXIS: xi, "thickness_mean_um": yi, "thickness_sd_um": si,
               "n_replicates": int(ni)} for xi, yi, si, ni in zip(x, y, sd, n_rep)]
    fits.append({"system": name, "code": code, "model": model.__name__, "param_1": popt[0],
                 "param_2": popt[1], "R2": r2, "x_definition": f"{X_AXIS} / {X_MAX}"})

xlabel = "Number of layers" if X_AXIS == "layers" else "Number of bilayers"
axs[0].set_ylabel("Thickness (µm)")
axs[2].set_ylabel("Thickness (µm)")
axs[2].set_xlabel(xlabel)
axs[3].set_xlabel(xlabel)

plt.tight_layout()
plt.savefig(FIGURE, dpi=600, bbox_inches="tight")
plt.close(fig)

pd.DataFrame(table).to_csv(OUT / "growth_data.csv", index=False)
pd.DataFrame(fits).to_csv(OUT / "growth_fits.csv", index=False)
print(pd.DataFrame(fits).round(3).to_string(index=False))
