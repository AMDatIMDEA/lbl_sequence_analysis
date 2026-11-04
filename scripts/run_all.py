"""
run_all.py
==========
Regenerates every file in processed_data/ and figures/ from raw_data/.

    cd scripts
    python run_all.py

Steps (each script can also be run on its own, in this order):
    01_thickness.py            raw thickness readings -> film thickness
    02_mechanical_properties.py  stress-strain curves + thickness -> film and architecture properties
    03_sequence_analysis.py    sequences -> composition, NCD, transitions, pairing bias, D_pair
    04_statistics.py           Welch ANOVA, Tukey, Pearson matrix, Figures S3-S4
    05_parent_variant.py       S -> NS / SS changes, Wilcoxon, Spearman vs NCD
    06_growth_test.py          growth of the four binary systems
    07_repeatability_test.py   10 Alg/Chi films
    08_pair_statistics.py      SI Section 9 (same as the notebook 08_pair_statistics.ipynb)
"""
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = ["01_thickness.py", "02_mechanical_properties.py", "03_sequence_analysis.py", "04_statistics.py",
         "05_parent_variant.py", "06_growth_test.py", "07_repeatability_test.py", "08_pair_statistics.py"]


def run(script, extra_env=None):
    env = {**os.environ, "MPLBACKEND": "Agg", **(extra_env or {})}
    print(f"\n=== {script} {extra_env or ''}")
    r = subprocess.run([sys.executable, script], cwd=HERE, env=env, capture_output=True, text=True)
    print(r.stdout[-1500:])
    if r.returncode:
        print(r.stderr[-3000:])
        sys.exit(f"{script} failed")


for s in STEPS:
    run(s)
print("\nAll done.")
