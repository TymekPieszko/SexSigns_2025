from pathlib import Path
import tskit
import numpy as np
import tqdm
import sys, pickle, json
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from functools import partial
import functions.phased as ph

# ----------------------------
# args
# ----------------------------
# python calc_delta_SINGER.py MP 0.1
# python calc_delta_SINGER.py CF 0.1
model = sys.argv[1]
bias = float(sys.argv[2])

# fixed params
rep_num = 100
window_size = 5000  # SINGER tree sequences are all 5kb

# directories
ts_dir = Path(f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/3.1.SINGER/")
checkpoint_file = Path(f"./delta/delta_{model}_bias_{bias}_SINGER.pkl")
out_file = Path(f"./delta/delta_{model}_bias_{bias}_SINGER.json")

# load or init checkpoint
if checkpoint_file.exists():
    with open(checkpoint_file, "rb") as f:
        total_reps = pickle.load(f)
else:
    total_reps = defaultdict(dict)

# ----------------------------
# worker
# ----------------------------
def calc_delta_master(rep_dir, window_size):
    """
    Compute δ for one replicate (averaged over MCMC samples inside rep_dir).
    Assumes each ts already contains exactly 2 individuals of length 5000.
    """
    mcmc_reps = []
    for ts_file in rep_dir.rglob("*.trees"):
        ts = tskit.load(str(ts_file))
        ts = ts.simplify(keep_input_roots=False)

        # genotypes + positions
        genos, pos = ph.get_biallelic_delta(ts)

        # tree-based runs only
        runs = ph.get_tree_runs(ts)

        delta = ph.run_delta(genos, pos, runs)
        mcmc_reps.append(delta)

    return np.nanmean(mcmc_reps)


# ----------------------------
# main loop
# ----------------------------
for dir in ts_dir.glob(f"SEX~*/REC~*/MUT~5e-07/BIAS~{bias}"):
    rec = float(dir.parts[-3].split("~")[1])
    sex = float(dir.parts[-4].split("~")[1])

    # skip if already have 100 reps
    if (
        sex in total_reps
        and rec in total_reps[sex]
        and len(total_reps[sex][rec]) == rep_num
    ):
        print(f"Skipping sex={sex}, rec={rec}")
        continue

    rep_dirs = sorted(dir.glob("*/ts"), key=lambda x: float(x.parts[-2]))[:rep_num]

    func = partial(calc_delta_master, window_size=window_size)

    with Pool(processes=cpu_count()) as pool:
        reps = list(
            tqdm.tqdm(
                pool.imap(func, rep_dirs),
                total=len(rep_dirs),
                desc=f"sex={sex}, rec={rec}"
            )
        )

    total_reps[f"{sex}_{rec}"] = reps

    with open(checkpoint_file, "wb") as f:
        pickle.dump(total_reps, f)

# ----------------------------
# save
# ----------------------------
print(total_reps)
with open(out_file, "w") as f:
    json.dump(total_reps, f, sort_keys=True)
