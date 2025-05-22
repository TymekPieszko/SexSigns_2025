from pathlib import Path
import sys, json
import numpy as np
import matplotlib.pyplot as plt
from functions_val import *

# python plot_validation.py MP Hi_windows_MP_mut_5e-07_inds_100
# python plot_validation.py CF Hi_windows_CF_mut_5e-07_inds_100
model = sys.argv[1]
in_name = sys.argv[2]
L = 1e6
TRACT = 5000
MUT = 5e-07
window_length = 10000

sim_data = f"/data/biol-bdelloids/scro4331/SexSigns_2025/validation/{in_name}.txt"
plot_file = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/validation/{in_name}.png"
)

with open(sim_data, "r") as f:
    sim_data = json.load(f)

sim_data = {
    SEX: {REC: np.nanmean(reps, axis=0) for REC, reps in REC_dict.items()}
    for SEX, REC_dict in sim_data.items()
}

SEX_lst = ["0.0"]
REC_lst = list(next(iter(sim_data.values())).keys())
combos = [(s, r) for s in SEX_lst for r in REC_lst]


fig, ax = plt.subplots(
    nrows=2, ncols=int(len(REC_lst) / 2), sharex=False, figsize=(15, 9)
)
ax = ax.flatten()
for i, combo in enumerate(combos):
    if model == "MP":
        positions = np.arange(0, L, window_length, dtype=int)
        hi_analytical = calc_hi(MUT, calc_gamma_MP(float(combo[1]), TRACT))
        hi_simulated = sim_data[str(combo[0])][str(combo[1])]
        ax[i].plot(
            positions,
            hi_simulated,
            color="cornflowerblue",
            marker="o",
            markersize=3,
            linestyle="",
        )
        ax[i].axhline(y=hi_analytical, color="black", linestyle="-")
        ax[i].set_title(r"$\gamma$ = " + f"{(float(combo[1]) * 5000):.2e}", fontsize=14)
        ax[i].set_ylabel(r"$H_{I}$", fontsize=12)
        ax[i].set_xlabel("Position (bp)", fontsize=12, labelpad=0.5)
    elif model == "CF":
        positions = np.arange(
            window_length, L + window_length, window_length, dtype=int
        ) - (window_length / 2)
        print(len(positions))
        hi_analytical = [
            calc_hi(MUT, calc_gamma_CF_NCI(float(combo[1]), pos)) for pos in positions
        ]
        print(len(hi_analytical))
        hi_simulated = sim_data[str(combo[0])][str(combo[1])]
        print(len(hi_simulated))
        ax[i].plot(
            positions,
            hi_simulated,
            color="cornflowerblue",
            marker="o",
            markersize=3,
            linestyle="",
        )
        ax[i].plot(positions, hi_analytical, color="black", linestyle="-")
        ax[i].set_xscale("log")
        ax[i].set_xlim(left=1e4, right=1e6)
        ax[i].set_title(
            r"$\gamma$ = "
            + f"{np.mean([calc_gamma_CF_NCI(float(combo[1]), pos) for pos in range(10**6)]):.2e}",
            fontsize=14,
        )
        ax[i].set_ylabel(r"$H_{I}$", fontsize=12)
        ax[i].set_xlabel("Position (bp)", fontsize=12, labelpad=0.5)

plt.tight_layout()
plt.savefig(plot_file)
plt.close()
