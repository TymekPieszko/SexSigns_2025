from pathlib import Path
import tskit, msprime, pyslim
import tqdm
from sexsigns_functions.calc import *
import sys
import numpy as np

# python calc_LD_curves.py MP 10
model = sys.argv[1]
NUM_INDS = int(sys.argv[2])
N_ANC = 1000
REC_ANC = 1.0e-6
MUT = 5e-7
ts_dir = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/0.ts/"
)

windows = np.arange(0, 100000, 1000)
for dir in ts_dir.glob("SEX~*/REC~*"):
    REC = float(dir.parts[-1].split("~")[1])
    SEX = float(dir.parts[-2].split("~")[1])
    # if SEX != 0.0:
    #     continue
    out_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/LD_curves/{model}/SEX_{SEX}_REC_{REC}.txt"
    reps = []
    for ts_file in tqdm.tqdm(dir.glob("*.trees")):
        ts = tskit.load(str(ts_file))
        ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC_ANC)
        ts = subsample_ts(ts, NUM_INDS)
        ts = msprime.sim_mutations(ts, rate=MUT)
        r2_windows = calc_r2_windows(ts, windows, 500)
        reps.append(r2_windows)
    reps_mean = np.nanmean(np.array(reps), axis=0)
    with open(out_file, "w") as f:
        for i in range(len(windows) - 1):
            f.write(f"{windows[i]} {reps_mean[i]}\n")
    f.close()
