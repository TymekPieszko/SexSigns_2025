import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import numpy as np
import itertools
import json, sys
from datetime import datetime
from sexsigns_functions.plot import calc_class_props, params
import functions.phased as ph
from functions.plot import combo_keys, sex_lst, rec_lst

# python plot_tcomp_combined_NEW.py MP
# python plot_tcomp_combined_NEW.py CF
model = sys.argv[1]
plot_file = Path(f"./heatmaps/tree_props_combined_{model}.png")
plotData_dir = f"./heatmaps/tree_props_plotData"
in_files = [
    f"{model}_sim_1.0",
    f"{model}_sim_0.1",
    # f"{model}_n_8_k_10_wind_5000_bias_0.1_clust",
    f"{model}_IQ_TREE_0.1",
    f"{model}_SINGER_0.1",
]
combos = list(itertools.product(in_files, range(3)))

if model == "MP":
    ylabels = ["5.0e-03", "5.0e-04", "5.0e-05", "5.0e-06", "0.0"]
elif model == "CF":
    ylabels = ["1.6e-01", "2.4e-02", "2.5e-03", "2.5e-04", "0.0"]
fig, ax = plt.subplots(
    3,
    4,
    sharex=True,
    figsize=(params["xfigsize"] * 4, params["yfigsize"] * 3),
)  # Keep size correspondence with the individual heatmaps
ax = ax.T.flatten()
for i, combo in enumerate(combos):
    in_file, j = combo
    print(in_file, j)
    sim_data = f"../stats/tcomp/{in_file}.json"
    with open(sim_data, "r") as f:
        sim_data = json.load(f)
    if "clust" in in_file:
        trees = ["CL", "AR", "SX"]
        tree_props = {k: np.mean([ph.calc_category_props(tree_dict)[trees[j]] for tree_dict in sim_data[k]]) for k in combo_keys}
        tree_props = {
            float(r): {
                float(s): tree_props.get(f"{s}_{r}", np.nan)
                for s in sex_lst
            }
            for r in rec_lst
        }
        tree_props = pd.DataFrame(tree_props).T.sort_index(ascending=False)
    else:
        trees = ["C", "A", "S"]
        tree_props = {
        float(SEX): {
            float(REC): np.nanmean(
                [calc_class_props(tree_dict=rep)[trees[j]] for rep in reps]
            )
            for REC, reps in REC_dict.items()
        }
        for SEX, REC_dict in sim_data.items()
        }
        tree_props = pd.DataFrame(tree_props).iloc[::-1]
        tree_props.sort_index(axis=0, inplace=True)
        tree_props.sort_index(axis=1, inplace=True)
        tree_props = tree_props.iloc[::-1, :-1]
    data_file = f"{plotData_dir}/{in_file}_{trees[j]}.csv"
    if not Path(data_file).exists():
        tree_props.to_csv(
            data_file,
            index=True,
            header=True,
            float_format="%.6f",
        )
    vmax = 1.0
    if trees[j] == "S" or trees[j] == "SX":
        vmax = 0.5
    sns.heatmap(tree_props, cmap="plasma", annot=False, vmin=0.0, vmax=vmax, ax=ax[i])
    ax[i].set_title(f"{in_file}", fontsize=params["title_font"], pad=16)
    ax[i].set_xticks(
        ticks=[0.5, 1.5, 3.5, 5.5, 7.5, 9.5],
        labels=["0.0", "1e-05", "1e-04", "1e-03", "1e-02", "1e-01"],
        fontsize=40,
        rotation=45,
    )
    ax[i].set_yticks(
        ticks=[0.5, 2.5, 4.5, 6.5, 7.5],
        labels=ylabels,
        fontsize=40,
        rotation=0,
    )
    colorbar = ax[i].collections[0].colorbar
    colorbar.ax.tick_params(labelsize=40)

for i in range(3):
    ax[i].set_ylabel(r"$\gamma$", fontsize=66, labelpad=16)
for i in [2, 5, 8, 11]:
    ax[i].set_xlabel(r"$\sigma$", fontsize=66, labelpad=16)
# plt.suptitle(
#     f"{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}", fontsize=params["title_font"]
# )
plt.tight_layout()
plt.subplots_adjust(hspace=0.18)
plt.savefig(plot_file, dpi=240)
