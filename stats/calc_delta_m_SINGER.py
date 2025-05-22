from pathlib import Path
import tskit
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import pickle, json
import numpy as np

# python calc_delta_m_SINGER.py MP 0.1
# python calc_delta_m_SINGER.py CF 0.1
model = sys.argv[1]
bias = float(sys.argv[2])
ts_dir = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/3.1.SINGER/"
)
checkpoint_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/delta_m/{model}_delta_m_bias_{bias}_SINGER.pkl"
out_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/delta_m/{model}_delta_m_bias_{bias}_SINGER.txt"

if Path(checkpoint_file).exists():
    with open(checkpoint_file, "rb") as f:
        total_reps = pickle.load(f)
else:
    total_reps = defaultdict(dict)
for dir in ts_dir.glob(f"SEX~*/REC~*/MUT~5e-07/BIAS~{bias}"):
    REC = float(dir.parts[-3].split("~")[1])
    SEX = float(dir.parts[-4].split("~")[1])
    if (
        SEX in total_reps
        and REC in total_reps[SEX]
        and len(total_reps[SEX][REC]) == 100
    ):
        print(f"Skipping {SEX} + {REC}")
        continue
    reps = []
    rep_dirs = sorted(dir.glob("*/ts"), key=lambda x: float(x.parts[-2]))
    for rep_dir in tqdm.tqdm(rep_dirs):
        # Note you are averageing across all MCMC samples across all windows; should not matter
        mcmc_reps = []
        for ts_file in rep_dir.rglob("*.trees"):
            ts = tskit.load(str(ts_file))
            ts = ts.simplify(
                keep_input_roots=False
            )  # Otherwise: ValueError: Cannot rank trees with unary nodes
            rep = calc_delta_m(ts)
            mcmc_reps.append(rep)
        mcmc_mean = np.nanmean(mcmc_reps)
        reps.append(mcmc_mean)
    # print(reps)
    total_reps[SEX][REC] = reps

    with open(checkpoint_file, "wb") as f:
        pickle.dump(total_reps, f)


print(total_reps)

with open(out_file, "w") as f:
    json.dump(total_reps, f, sort_keys=True)
