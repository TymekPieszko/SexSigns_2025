import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import json
import sys
from sexsigns_functions.plot import calc_aic, params


# python plot_Fis_analysis.py MP
model = sys.argv[1]
obs = 0.0
# MUT_lst = [5e-09, 5e-08, 5e-07]
INDS_lst = [5, 10, 15]
# combos = [(i, j) for i in MUT_lst for j in INDS_lst]
plot_file = f"./Fis_analysis/{model}_obs_{obs}.png"


def calc_likelihood(sim_data, obs):
    likelihoods = pd.DataFrame(
        columns=list(sim_data.keys()),
        index=list(sim_data[list(sim_data.keys())[0]].keys())[::-1],
    )
    for i, REC_RATE in enumerate(likelihoods.index):
        for j, SEX_FREQ in enumerate(likelihoods.columns):
            mean = np.mean(sim_data[str(SEX_FREQ)][str(REC_RATE)])
            std = np.std(sim_data[str(SEX_FREQ)][str(REC_RATE)])
            likelihood = stats.norm.pdf(obs, loc=mean, scale=std)
            if likelihood == 0:
                log_likelihood = -np.inf
            else:
                log_likelihood = np.log(likelihood)
            likelihoods.iloc[i, j] = log_likelihood
    return likelihoods


nrows = 1
ncols = 3
fig, ax = plt.subplots(
    nrows,
    ncols,
    sharex=True,
    sharey=True,
    figsize=(30, 10),
)
ax = ax.T.flatten()
for i, INDS in enumerate(INDS_lst):
    in_name = f"Fis_{model}_mut_5e-07_inds_{INDS}"
    data_path = f"../stats/Fis/{in_name}.txt"
    with open(data_path, "r") as f:
        sim_data = json.load(f)
    # print(sim_data)
    likelihoods = calc_likelihood(sim_data, obs).astype(float)[::-1]
    aic = likelihoods.transform(calc_aic, n_params=3)
    aic = aic.to_numpy()
    delta_aic = aic - np.min(aic)
    rel_likelihoods = np.exp(-0.5 * delta_aic)
    aic_weights = rel_likelihoods / np.sum(rel_likelihoods)
    # print(aic_weights)
    # Plot the Akaike weights
    ax[i].plot(
        likelihoods.index,
        aic_weights,
        marker="o",
        color="black",
        linewidth=4,
        markersize=8,
    )
    top_indices = np.argsort(aic_weights)[-3:]
    # for x_val, y_val in zip(likelihoods.index, aic_weights):
    #     ax[i].annotate(
    #         f"{y_val}",
    #         xy=(x_val, y_val),
    #         xytext=(0, 5),
    #         textcoords="offset points",
    #         ha="center",
    #         fontsize=12,
    #     )
    ax[i].set_title(f"n = {str(INDS)}", fontsize=28)
    if model == "MP":
        labels = [
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
        labels = [
            "0.0",
            "2.5e-04",
            "7.9e-04",
            "2.5e-03",
            "7.7e-03",
            "2.4e-02",
            "6.8e-02",
            "1.6e-01",
        ]
    ax[i].set_xticks(
        ticks=np.arange(len(likelihoods.index)),
        labels=labels,
        fontsize=params["tick_font"],
        rotation=45,
    )
    ax[i].set_ylim(0, 1)
    ax[i].tick_params(axis="y", labelsize=params["tick_font"])

    aic_weights = sorted(aic_weights, reverse=True)
    print(in_name)
    print("Top 3:", aic_weights[0:3])
    # Find the x range covering 90% of the weight
    # pairs = sorted(
    #     zip(likelihoods.index, aic_weights), key=lambda x: x[1], reverse=True
    # )
    # covered = []
    # total = 0
    # for x_val, w in pairs:
    #     covered.append(x_val)
    #     total += w
    #     if total >= 0.9:
    #         break
    # print(f"90% weight range: {min(covered)} to {max(covered)}")
ax[0].set_ylabel(r"$A_{\mathit{w}}$", fontsize=params["label_font"], labelpad=16)
for i in [0, 1, 2]:
    ax[i].set_xlabel(r"$\gamma$", fontsize=params["label_font"], labelpad=16)
plt.tight_layout()
plt.subplots_adjust(hspace=0.1, wspace=0.05)
plt.savefig(plot_file)
