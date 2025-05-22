import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import numpy as np
import json, sys
from datetime import datetime
from sexsigns_functions.plot import params

# python plot_95_limits.py MP 0.001
model = sys.argv[1]
target_sex = sys.argv[2]
if model == "MP":
    rec_values = [
        "0.0",
        "5e-06",
        "1.6e-05",
        "5e-05",
        "1.6e-04",
        "5e-04",
        "1.6e-03",
        "5e-03",
    ]
elif model == "CF":
    rec_values = [
        "0.0",
        "2.5e-04",
        "7.9e-04",
        "2.5e-03",
        "7.7e-03",
        "2.4e-02",
        "6.8e-02",
        "1.6e-01",
    ]
stats = [
    "Hi",
    "Fis",
    "LD_0.1000",
    "S_1_0",
    "S_0_1",
    "S_IQ_TREE",
    "S_SINGER",
    "D_1_0",
    "D_0_1",
    "D_IQ_TREE",
    "D_SINGER",
]
stats_pretty = [
    r"$H_I$",
    r"$F_{IS}$",
    r"$R^2$" + " < 1000 bp",
    r"$p(SX)$",
    r"$p(SX)$" + " + 0.1",
    r"$p(SX)$" + " + IQ_TREE",
    r"$p(SX)$" + " + SINGER",
    r"$\Delta_m$",
    r"$\Delta_m$" + " + 0.1",
    r"$\Delta_m$" + " + IQ_TREE",
    r"$\Delta_m$" + " + SINGER",
]
in_dir = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/variation_analysis/95_limits_SEX_{target_sex}_{model}/"
)
plot_file = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/variation_analysis/plots/95_limits_SEX_{target_sex}_{model}.png"
)
########################################
colors = plt.cm.viridis(np.linspace(0, 1, 8))
print(colors)
fig, ax = plt.subplots(
    3,
    4,
    sharex=True,
    sharey=True,
    figsize=(params["xfigsize"] * 4, params["yfigsize"] * 3),
)  # Keep size correspondence with the individual heatmaps
ax = ax.flatten()
for i, stat in enumerate(stats):
    i2 = i
    if i > 2:
        i2 = i + 1
    file = in_dir / f"{stat}.txt"
    data = pd.read_csv(file, sep="\t", header=0, index_col=0)
    # print(data)
    for j in range(len(data.index)):
        rec_rate_data = data.iloc[j, :]
        ax[i2].plot(
            list(range(len(rec_rate_data))), rec_rate_data, color=colors[j], lw=6
        )
        ax[i2].set_ylim(0, 1)
        ax[i2].set_title(
            f"{stats_pretty[i]}",
            fontsize=60,
            pad=33,
        )
        ax[i2].set_xticks(
            ticks=[0, 1, 3, 5, 7, 9],
            labels=["0.0", "1e-05", "1e-04", "1e-03", "1e-02", "1e-01"],
            fontsize=42,
            rotation=45,
        )
        ax[i2].set_yticks(
            ticks=np.arange(0, 1.1, 0.2),
            labels=["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"],
            fontsize=42,
            rotation=0,
        )
        ax[i2].axvspan(-0.5, 0.5, color="gray", alpha=0.05)
        ax[i2].axvspan(0.5, 5.5, color="pink", alpha=0.05)
        ax[i2].axvspan(5.5, 10.5, color="crimson", alpha=0.05)
        ax[i2].axvline(
            x=5,
            color="gray",
            lw=6,
            ls="--",
        )

ax[3].axis("off")
legend_lines = [
    plt.Line2D([0], [0], color=colors[j], lw=6) for j in range(len(rec_values))
]

# Add the legend to ax[3]
ax[3].legend(
    legend_lines,
    rec_values,
    title=r"$\gamma$",
    fontsize=50,
    title_fontsize=52,
    loc="center",
    frameon=True,
    ncol=2,
)
fig.text(
    0.01,
    0.5,
    "Prop. within 95% interval (" + r"$\sigma$" + f" = {float(target_sex)})",
    va="center",
    rotation="vertical",
    fontsize=68,
)
fig.text(0.5, 0.01, r"$\sigma$", ha="center", fontsize=68)
plt.tight_layout(rect=(0.04, 0.04, 1, 1))
plt.subplots_adjust(hspace=0.22, wspace=0.06)
plt.savefig(plot_file)
