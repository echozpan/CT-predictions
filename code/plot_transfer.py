# Plot cross-country transfer penalties and calibration by CT quartile.
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
res = Path(sys.argv[1])
fig_dir = res / "figures"
fig_dir.mkdir(exist_ok=True)
order = ["DEU", "DNK", "FIN", "FRA", "KOR", "LUX", "PRT", "USA"]
names = {"DEU": "Germany", "DNK": "Denmark", "FIN": "Finland", "FRA": "France", "KOR": "Korea", "LUX": "Luxembourg", "PRT": "Portugal", "USA": "United States"}

                                                                         
pen = pd.read_csv(res / "transfer_penalties.csv")
pooled = pd.read_csv(res / "transfer_matrix_pooled.csv")
for fs, tag in [("All predictors", "all"), ("Excluding PV1CIL", "nocil")]:
    sub = pen[(pen.model == "Lasso") & (pen.feature_set == fs) & (pen.source.isin(order))]
    mat = sub.pivot(index="source", columns="target", values="relative_penalty_percent").reindex(index=order, columns=order)
    rm = pooled[(pooled.model == "Lasso") & (pooled.feature_set == fs) & (pooled.is_within) & (pooled.source.isin(order))].set_index("target").RMSE
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    vmax = np.nanmax(mat.values)
    im = ax.imshow(mat.values, cmap="Greys", vmin=0, vmax=vmax)
    ax.set_xticks(range(8)); ax.set_xticklabels([names[c] for c in order], rotation=40, ha="right")
    ax.set_yticks(range(8)); ax.set_yticklabels([names[c] for c in order])
    ax.set_xlabel("Target sample (tested on)"); ax.set_ylabel("Source sample (trained on)")
    for i, s in enumerate(order):
        for j, t in enumerate(order):
            if s == t:
                ax.text(j, i, f"{rm[t]:.0f}", ha="center", va="center", fontsize=8, color="black", fontweight="bold")
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, ec="black", lw=1.2))
            else:
                v = mat.values[i, j]
                ax.text(j, i, f"{v:+.0f}%", ha="center", va="center", fontsize=8, color="white" if v > 0.6 * vmax else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("RMSE increase over local benchmark (%)")
    ax.set_title("Cross-country transfer penalty (lasso" + ("" if tag == "all" else ", no CIL") + ")")
    fig.tight_layout(); fig.savefig(fig_dir / f"transfer_matrix_{tag}.png", dpi=300); plt.close(fig)

                                                                                                              
fs_file = res / "fairness_subgroup_errors.csv"
if fs_file.exists():
    f = pd.read_csv(fs_file)
    f = f[(f.dimension == "Observed CT quartile") & (f.feature_set == "All predictors") & (f.model == "Lasso")]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.1), sharey=True)
    for ax, code in zip(axes, ["KOR", "DNK"]):
        for ctx, ls, mk in [("Within-country", "-", "o"), ("Europe-trained external", "--", "s")]:
            s = f[(f.target_code == code) & (f.context == ctx)].sort_values("subgroup")
            ax.plot(range(1, 5), s.mean_error, ls, marker=mk, color="black", label=ctx.replace("Europe-trained external", "Europe-trained"), markerfacecolor="white" if mk == "s" else "black")
        ax.axhline(0, color="gray", lw=.8)
        ax.set_xticks(range(1, 5)); ax.set_xticklabels(["Q1\nlowest", "Q2", "Q3", "Q4\nhighest"])
        ax.set_title({"KOR": "Korea", "DNK": "Denmark"}[code]); ax.set_xlabel("Observed CT quartile")
    axes[0].set_ylabel("Mean error (predicted − observed)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(fig_dir / "calibration_by_quartile.png", dpi=300); plt.close(fig)
print("figures written to", fig_dir)
