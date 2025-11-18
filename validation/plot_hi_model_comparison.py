from pathlib import Path
from functions_val import *
import numpy as np
import matplotlib.pyplot as plt

L = 1e6
TRACT = 5000
MUT = 5e-07
n_lst = [
    0.001,
    0.01,
    0.1,
    1,
    10,
]  # Mean of the Poisson-distributed number of crossovers

plot_file = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/validation/model_comparison.png"
)

fig, ax = plt.subplots(figsize=(10, 6))
plt.axhline(y=1, color="red", linewidth=2)
for i, n in enumerate(n_lst):
    hi_CF_CI = calc_hi(MUT, calc_gamma_CF_CI(REC=n, pos=1))
    hi_CF_NCI = calc_hi(MUT, calc_gamma_CF_NCI(REC=n, pos=1))
    ratio = hi_CF_CI / hi_CF_NCI
    ax.hlines(y=ratio, xmin=i - 0.3, xmax=i + 0.3, color="black", linewidth=3)
plt.xticks(range(len(n_lst)), n_lst)
plt.xlabel(r"$\bar{n}$", fontsize=12)
plt.ylabel(r"$\dfrac{H_{I} \mid CI}{H_{I} \mid NCI}$", fontsize=12)
plt.ylim(0.6, 1.4)
plt.tight_layout()
plt.savefig(plot_file)
