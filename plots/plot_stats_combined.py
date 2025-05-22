import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime
import json, sys
from sexsigns_functions.plot import params

# python plot_stats_combined_per_model.py MP
# python plot_stats_combined_per_model.py CF
model = sys.argv[1]
in_files = [
    f"Hi_{model}_mut_5e-07_inds_100",
    f"Fis_{model}_mut_5e-07_inds_100",
    f"LD_0-1000_{model}_mut_5e-07_inds_100",
    # f"LD_9000-10000_{model}_mut_5e-07_inds_100",
]
plot_file = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/plots/heatmaps/simulated/stats_combined_{model}.png"
)
nrows = 1
ncols = len(in_files)
fig, ax = plt.subplots(
    nrows,
    ncols,
    sharex=True,
    figsize=(36, 10),
)
ax = ax.T.flatten()
for i, in_file in enumerate(in_files):
    sim_data = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/{in_file.split(f'_{model}')[0]}/{in_file}.txt"
    print(sim_data)
    with open(sim_data, "r") as f:
        sim_data = json.load(f)
    means = {
        key1: {key2: np.mean(val2) for key2, val2 in val1.items()}
        for key1, val1 in sim_data.items()
    }
    means = pd.DataFrame(means).iloc[::-1]
    means = means.iloc[:, :-1]
    # print(means)
    sns.heatmap(means, cmap="plasma", annot=False, ax=ax[i])
    ax[i].set_title(f"{in_file}", fontsize=params["title_font"], pad=16)
    ax[i].set_xticks(
        ticks=[0.5, 1.5, 3.5, 5.5, 7.5, 9.5],
        labels=["0.0", "1.0e-05", "1.0e-04", "1.0e-03", "1.0e-02", "1.0e-01"],
        fontsize=params["tick_font"],
        rotation=45,
    )
    colorbar = ax[i].collections[0].colorbar
    colorbar.ax.tick_params(labelsize=params["tick_font"])


if model == "MP":
    ylabels = ["5.0e-03", "5.0e-04", "5.0e-05", "5.0e-06", "0.0"]
elif model == "CF":
    ylabels = ["1.6e-01", "2.4e-02", "2.5e-03", "2.5e-04", "0.0"]
for i in range(ncols):
    ax[i].set_yticks(
        ticks=[0.5, 2.5, 4.5, 6.5, 7.5],
        labels=ylabels,
        fontsize=params["tick_font"],
        rotation=0,
    )

for i in range(nrows):
    ax[i].set_ylabel(r"$\gamma$", fontsize=params["label_font"], labelpad=16)
for i in range(ncols):
    ax[i].set_xlabel(r"$\mathit{\sigma}$", fontsize=params["label_font"], labelpad=16)
plt.suptitle(
    f"{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}", fontsize=params["title_font"]
)
plt.tight_layout()
plt.subplots_adjust(hspace=0.18)
plt.savefig(plot_file)
