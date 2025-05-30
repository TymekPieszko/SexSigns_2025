from pathlib import Path
import pandas as pd
import numpy as np
import json
import sys
from sexsigns_functions.plot import calc_class_props
import matplotlib.pyplot as plt

# python plot_tree_lines.py MP MP_sim_1.0 MP_sim_0.1 S
# python plot_tree_lines.py MP MP_sim_0.1 MP_IQ_TREE_0.1 S
# python plot_tree_lines.py MP MP_sim_0.1 MP_SINGER_0.1 S
model = sys.argv[1]
tree = sys.argv[4]
if tree == "C":
    tree_index = 0
elif tree == "A":
    tree_index = 1
elif tree == "S":
    tree_index = 2
if model == "MP":
    gammas = [
        "0.0",
        "5e-06",
        "1.6e-05",
        "5e-05",
        "1.6e-04",
        "5e-04",
        "1.6e-03",
        "5e-03",
    ][::-1]
elif model == "CF":
    gammas = [
        "0.0",
        "2.5e-04",
        "7.9e-04",
        "2.5e-03",
        "7.7e-03",
        "2.4e-02",
        "6.8e-02",
        "1.6e-01",
    ][::-1]
plot_file = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/plots/tree_prop_lines/{sys.argv[2]}_{sys.argv[3]}_{tree}_lines.png"
)
REC_lst = [0.0, 1.0e-09, 3.162e-09, 1.0e-08, 3.162e-08, 1.0e-07, 3.162e-07, 1.0e-06][
    ::-1
]
SEX_lst = [
    0.0,
    1.0e-05,
    3.162e-05,
    1.0e-04,
    3.162e-04,
    1.0e-03,
    3.162e-03,
    1.0e-02,
    3.162e-02,
    0.1,
    3.162e-01,
]
xlabels = [
    "0.0",
    "1.0e-05",
    "3.162e-05",
    "1.0e-04",
    "3.162e-04",
    "1.0e-03",
    "3.162e-03",
    "1.0e-02",
    "3.162e-02",
    "0.1",
    "3.162e-01",
]
combos = [(s, r) for s in SEX_lst for r in REC_lst]


### PLOT ####
fig, ax = plt.subplots(nrows=len(REC_lst), ncols=1, sharex=True, figsize=(16, 22))
ax = ax.T.flatten()
line_data = {i: [None, None] for i in range(len(REC_lst))}
# print(line_data)
for i in [2, 3]:
    in_name = sys.argv[i]
    sim_data = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/tree_composition/{in_name}.txt"

    with open(sim_data, "r") as f:
        sim_data = json.load(f)

    class_props = {
        SEX: {
            REC: np.nanmean(
                [list(calc_class_props(tree_dict=rep).values()) for rep in reps], axis=0
            )
            for REC, reps in REC_dict.items()
        }
        for SEX, REC_dict in sim_data.items()
    }

    for j, REC in enumerate(REC_lst):
        CAS_values = np.array([class_props[str(SEX)][str(REC)] for SEX in SEX_lst])
        ax[j].plot(
            CAS_values[:, tree_index],
            color="black",
            linestyle="-" if i == 2 else "--",
            linewidth=2,
        )
        ax[j].set_title(r"$\gamma$" + f" = {gammas[j]}", fontsize=20)
        ax[j].set_ylim(0, 1)
        ax[j].set_xticks(
            np.arange(len(SEX_lst)),
            labels=xlabels,
            fontsize=12,
            rotation=45,
        )
        line_data[j][i - 2] = CAS_values[:, 2]

        # if i == 3:
        #     mask = line_data[j][1] > line_data[j][0]
        #     ax[j].fill_between(
        #         range(len(line_data[j][0])),
        #         0,
        #         1,
        #         where=mask,
        #         color="lightblue",
        #         alpha=0.5,
        #     )
fig.text(
    0.015,
    0.5,
    "p(SX)",
    ha="center",
    rotation="vertical",
    fontsize=26,
)
fig.text(
    0.5,
    0.015,
    r"$\sigma$",
    ha="center",
    fontsize=26,
)
plt.tight_layout(rect=(0.03, 0.03, 1, 1))
plt.savefig(plot_file, dpi=320)
