from pathlib import Path
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys
from datetime import datetime
from sexsigns_functions.plot import calc_class_props, plot_heatmap, params

# python plot_LD_curves.py MP

model = sys.argv[1]

dir = f"../stats/LD_curves/{model}/"
plot_file = f"./LD_curves/{model}.png"

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
    1.0,
]
REC_lst = [0.0, 1.0e-09, 3.162e-09, 1.0e-08, 3.162e-08, 1.0e-07, 3.162e-07, 1.0e-06][
    ::-1
]
GAMMA_MP_lst = [
    "0.0",
    "5e-06",
    "1.6e-05",
    "5e-05",
    "1.6e-04",
    "5e-04",
    "1.6e-03",
    "5e-03",
][::-1]
GAMMA_CF_lst = [
    "0.0",
    "2.5e-04",
    "7.9e-04",
    "2.5e-03",
    "7.7e-03",
    "2.4e-02",
    "6.8e-02",
    "1.6e-01",
][::-1]
gamma_dict = {rec: (GAMMA_MP_lst[i], GAMMA_CF_lst[i]) for i, rec in enumerate(REC_lst)}
print(gamma_dict)
combos = [(s, r) for s in SEX_lst for r in REC_lst]

fig, ax = plt.subplots(
    nrows=len(REC_lst), ncols=len(SEX_lst), sharex=True, figsize=(28, 20)
)
ax = ax.T.flatten()
for i, combo in enumerate(combos):
    SEX = combo[0]
    REC = combo[1]
    r2_windows = np.loadtxt(f"{dir}/SEX_{SEX}_REC_{REC}.txt")
    ax[i].plot(r2_windows[:, 0], r2_windows[:, 1], color="black")
    if model == "MP":
        j = 0
    elif model == "CF":
        j = 1
    ax[i].set_title(f"$\\sigma$={SEX}, $\\gamma$={gamma_dict[REC][j]}")
    ax[i].set_xticks(np.arange(0, 100001, 10000))
    ax[i].set_xticklabels(["0", "", "", "", "", "", "", "", "", "", "100 kb"])
    # group_size = 5
    # x_means = []
    # y_means = []
    # for start in range(0, len(r2_windows), group_size):
    #     subset = r2_windows[start : start + group_size]
    #     if len(subset) < group_size:
    #         break
    #     x_means.append(subset[:, 0].mean())
    #     y_means.append(subset[:, 1].mean())
    # ax[i].plot(x_means, y_means, color="black", lw=2)
# for i in range(len(REC_lst)):
#     ax[i].set_ylabel(
#         r"R²",
#         fontsize=params["label_font"],
#         labelpad=16,
#     )
# for i in [7, 15, 23, 31, 39, 47]:
#     ax[i].set_xlabel(
#         "Distance (bp)",
#         fontsize=params["label_font"],
#         labelpad=16,
#     )
fig.text(
    0.01,
    0.5,
    r"$\mathit{R^2}$",
    va="center",
    rotation="vertical",
    fontsize=params["label_font"],
)
fig.text(0.5, 0.01, "Distance (bp)", ha="center", fontsize=params["label_font"])
# for i in [1, 3, 5]:
#     ax[i].set_xlabel(r"$\mathit{\sigma}$", fontsize=params["label_font"], labelpad=16)
plt.tight_layout(rect=(0.04, 0.04, 1, 1))
plt.savefig(plot_file, dpi=320)
plt.close()
