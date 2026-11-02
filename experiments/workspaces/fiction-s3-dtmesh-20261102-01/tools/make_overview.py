#!/usr/bin/env python3
"""Generate the two-panel S2/S3 overview figure for the blog.

Left: S2 spatial refinement (M1-M6).
Right: S3 time-step difference for M2/M3/M4 (from analysis-summary.json).

Usage:
  /home/steven/miniconda3/envs/elastogen/bin/python tools/make_overview.py
"""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]  # .../fiction
FIG = ROOT / "figures/validation/s2-s3-refinement-20261102-01"
S2 = ROOT / "figures/validation/s2-s3-refinement-20261102-01/s2-refinement.csv"
S3 = FIG / "analysis-summary.json"

s2 = list(csv.DictReader(S2.open()))
s3 = json.load(S3.open())
errs = s3["time_errors"]

fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.8))
ax = axes[0]
x = np.array([float(r["nx"]) for r in s2])
y = np.array([float(r["error_percent"]) for r in s2])
ax.loglog(x, y, "o-", color="#1f77b4", lw=1.5, ms=5)
for xi, yi, r in zip(x, y, s2):
    ax.annotate(r["mesh"], (xi, yi), textcoords="offset points", xytext=(4, 4), fontsize=8)
ax.set_xlabel("nx (elements along the bar)")
ax.set_ylabel("tip error vs beam theory [%]")
ax.set_title("S2: spatial refinement")
ax.grid(True, which="both", alpha=0.25)

ax = axes[1]
colors = {"M2": "#1f77b4", "M3": "#d62728", "M4": "#2ca02c"}
for r in errs:
    xs = np.array([0.005, 0.0025])
    ys = np.array([r["e_max_0.005_vs_0.0025_m"], r["e_max_0.0025_vs_0.00125_m"]]) * 1000
    ax.plot(xs, ys, "o-", color=colors.get(r["mesh"], "gray"), lw=1.5, ms=5,
            label=f"{r['mesh']} (p={r['observed_order_p']:.2f})")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([0.005, 0.0025])
ax.set_xticklabels(["0.005", "0.0025"])
ax.set_xlabel("dt [s]")
ax.set_ylabel("max tip difference to next finer dt [mm]")
ax.set_title("S3: time-step refinement")
ax.grid(True, which="both", alpha=0.25)
ax.legend(frameon=False, fontsize=8)

fig.tight_layout()
fig.savefig(FIG / "s2-s3-overview.png", dpi=200)
fig.savefig(FIG / "s2-s3-overview.svg")
print("wrote", FIG / "s2-s3-overview.png")
