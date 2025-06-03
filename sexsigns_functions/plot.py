# from scipy.interpolate import RegularGridInterpolator
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

# import itertools


# def interpolate(df, inter_num):
#     REC_RATE_lst = list(df.index.astype(float))
#     SEX_FREQ_lst = list(df.columns.astype(float))
#     interp = RegularGridInterpolator((REC_RATE_lst, SEX_FREQ_lst), df.values)

#     # Generate intermediate rec rates
#     REC_RATE_lst[-1] = 1.0e-10
#     REC_RATE_lst = list(
#         np.concatenate(
#             [
#                 np.geomspace(r1, r2, inter_num, endpoint=False)
#                 for r1, r2 in zip(REC_RATE_lst[:-1], REC_RATE_lst[1:])
#             ]
#         )
#     )
#     REC_RATE_lst.append(0.0)

#     # Generate intermediate sex frequencies
#     SEX_FREQ_lst[0] = 1e-06
#     SEX_FREQ_lst = list(
#         np.concatenate(
#             [
#                 np.geomspace(r1, r2, inter_num, endpoint=False)
#                 for r1, r2 in zip(SEX_FREQ_lst[:-1], SEX_FREQ_lst[1:])
#             ]
#         )
#     )
#     SEX_FREQ_lst.append(1.0)

#     coordinates = list(itertools.product(REC_RATE_lst, SEX_FREQ_lst))
#     interp_values = interp(coordinates)
#     interp_values = interp_values.reshape((len(REC_RATE_lst), len(SEX_FREQ_lst)))

#     return interp_values


def calc_likelihood(sim_data, obs):
    likelihoods = pd.DataFrame(
        columns=list(sim_data.keys())[:-1],
        index=list(sim_data[list(sim_data.keys())[0]].keys())[::-1],
    )
    print(likelihoods)
    for i, REC in enumerate(likelihoods.index):
        for j, SEX in enumerate(likelihoods.columns):
            print(REC, SEX_RATE)
            mean = np.mean(sim_data[str(SEX_FREQ)][str(REC_RATE)])
            std = np.std(sim_data[str(SEX_FREQ)][str(REC_RATE)])
            likelihood = stats.norm.pdf(obs, loc=mean, scale=std)
            if likelihood == 0:
                log_likelihood = -np.inf
            else:
                log_likelihood = np.log(likelihood)
            likelihoods.iloc[i, j] = log_likelihood
    return likelihoods


def calc_aic(log_likelihood, n_params):
    return 2 * n_params - 2 * log_likelihood


# calc_posterior()


def calc_class_props(tree_dict):
    C = (
        sum(tree_dict["C"].values())
        + tree_dict["S"]["b0213_02"]
        + tree_dict["S"]["b0213_13"]
        + tree_dict["S"]["b0312_03"]
        + tree_dict["S"]["b0312_12"]
    )
    S = (
        sum(tree_dict["S"].values())
        - tree_dict["S"]["b0213_02"]
        - tree_dict["S"]["b0213_13"]
        - tree_dict["S"]["b0312_03"]
        - tree_dict["S"]["b0312_12"]
        - tree_dict["S"]["b0123_eq"]
    )
    A = sum(tree_dict["A"].values()) + tree_dict["S"]["b0123_eq"]
    total = C + A + S
    if total == 0:
        return {"C": np.nan, "A": np.nan, "S": np.nan}
    return {"C": C / total, "A": A / total, "S": S / total}


def plot_heatmap(
    df,
    figsize,
    plot_file,
    title,
    xlab,
    ylab,
    label_font,
    tick_font,
    title_font,
    vmin,
    vmax,
    annot,
):
    plt.figure(figsize=figsize)
    ax = sns.heatmap(df, cmap="plasma", vmin=vmin, vmax=vmax, annot=annot)
    plt.ylabel(ylab, fontsize=label_font, labelpad=16)
    plt.xlabel(xlab, fontsize=label_font, labelpad=16)
    plt.xticks(
        ticks=[0.5, 1.5, 3.5, 5.5, 7.5, 9.5],
        labels=["0.0", "1e-05", "1e-04", "1e-03", "1e-02", "1e-01"],
        fontsize=tick_font,
        rotation=45,
    )
    plt.yticks(
        ticks=[0.5, 2.5, 4.5, 6.5, 7.5],
        labels=["0.0", "1e-09", "1e-08", "1e-07", "1e-06"][::-1],
        fontsize=tick_font,
        rotation=0,
    )
    plt.title(f"{title}", fontsize=title_font, pad=16)
    colorbar = ax.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=tick_font)
    plt.tight_layout()
    plt.savefig(plot_file)


params = {
    "xfigsize": 14,
    "yfigsize": 10,
    "label_font": 44,
    "tick_font": 26,
    "title_font": 20,
}

# See if there is a reason to add this here!!!
# def calc_gamma(pos, REC):
#     return (1 / 3) * (1 - np.exp(-1.5 * pos * REC))
