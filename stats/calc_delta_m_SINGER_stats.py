from pathlib import Path
import tskit, msprime, pyslim
import tqdm
from calc_functions import *
import sys
from collections import defaultdict
import json
import numpy as np

# python calc_delta_m_SINGER_stats.py MP 0.1
model = sys.argv[1]
bias = float(sys.argv[2])
ts_dir = Path(
    f"../sim_pipeline/sim_output/{model}/3.1.SINGER/"
)
out_file = f"./delta_m/{model}_delta_m_bias_{bias}_SINGER_stats.txt"

total_reps = defaultdict(dict)
for dir in ts_dir.glob(f"SEX~*/REC~*/MUT~5e-07/BIAS~{bias}"):
    REC = float(dir.parts[-3].split("~")[1])
    SEX = float(dir.parts[-4].split("~")[1])
    reps = []
    for rep_dir in tqdm.tqdm(dir.glob("*/ts")):
        target_trees = 0
        for ts_file in rep_dir.rglob("*.trees"):
            ts = tskit.load(str(ts_file))
            ts = ts.simplify(
                keep_input_roots=False
            )  # Otherwise: ValueError: Cannot rank trees with unary nodes
            for tree in ts.trees():
                rank = tree.rank()
                if rank in [(4, 1), (4, 2)]:
                    target_trees += 1
        reps.append(target_trees)
    print(SEX, REC, np.mean(reps))
    total_reps[SEX][REC] = np.mean(reps)

print(total_reps)

with open(out_file, "w") as f:
    for SEX in total_reps.keys():
        for REC in total_reps[SEX].keys():
            f.write(f"{SEX} {REC} {total_reps[SEX][REC]}\n")
