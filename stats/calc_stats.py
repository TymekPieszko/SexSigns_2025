from pathlib import Path
import tskit, msprime, pyslim
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import json

# python calc_stats.py MP Hi 100
# python calc_stats.py MP Fis 100
# python calc_stats.py MP LD_0-1000 100
model = sys.argv[1]
stat = sys.argv[2]
NUM_INDS = int(sys.argv[3])
N_ANC = 1000
REC_ANC = 1.0e-6
MUT = 5e-07
ts_dir = Path(
    f"../sim_pipeline/sim_output/{model}/0.ts/"
)
out_file = f"./{stat}/{stat}_{model}_mut_{MUT}_inds_{NUM_INDS}.txt"

total_reps = defaultdict(dict)
for dir in ts_dir.glob("SEX~*/REC~*"):
    REC = float(dir.parts[-1].split("~")[1])
    SEX = float(dir.parts[-2].split("~")[1])
    reps = []
    ts_files = sorted(dir.glob("*.trees"), key=lambda x: float(x.stem))
    for ts_file in tqdm.tqdm(ts_files):
        ts = tskit.load(str(ts_file))
        ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC_ANC)
        ts = subsample_ts(ts, NUM_INDS)
        ts = msprime.sim_mutations(ts, rate=MUT)
        print(ts.num_sites)
        if stat == "Hi":
            rep = calc_hi(ts)
        elif stat == "Fis":
            rep = calc_fis(ts)
        elif stat == "LD_0-1000":
            rep = calc_r2(ts, 0, 1000)
        elif stat == "LD_9000-10000":
            rep = calc_r2(ts, 9000, 10000)
        reps.append(rep)
    total_reps[SEX][REC] = reps

print(total_reps)
with open(out_file, "w") as f:
    json.dump(total_reps, f, sort_keys=True)
