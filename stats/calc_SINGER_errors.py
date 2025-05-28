from pathlib import Path
import tqdm
from collections import defaultdict
import sys
import numpy as np

# python calc_SINGER_errors.py MP 0.1
model = sys.argv[1]
bias = float(sys.argv[2])
ts_dir = Path(
    f"../sim_pipeline/sim_output/{model}/3.1.SINGER/"
)
out_file = f"./SINGER_errors/{model}_errors.txt"

total_reps = defaultdict(dict)
for dir in ts_dir.glob(f"SEX~*/REC~*/MUT~5e-07/BIAS~{bias}"):
    REC = float(dir.parts[-3].split("~")[1])
    SEX = float(dir.parts[-4].split("~")[1])
    print(SEX, REC)
    reps = []
    for rep_dir in tqdm.tqdm(dir.glob("*/ts")):
        mcmc_reps = []
        mcmc_num = len(list(rep_dir.rglob("*.trees")))
        reps.append(mcmc_num)
    total_reps[SEX][REC] = np.nanmean(reps)

print(total_reps)

with open(out_file, "w") as f:
    for SEX in total_reps.keys():
        for REC in total_reps[SEX].keys():
            f.write(f"{SEX}\t{REC}\t{total_reps[SEX][REC] / 200}\n")
