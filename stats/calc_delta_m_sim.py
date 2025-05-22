from pathlib import Path
import tskit, msprime, pyslim
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import json

# python calc_delta_m_sim.py CF 1.0
# python calc_delta_m_sim.py MP 1.0
model = sys.argv[1]
bias = float(sys.argv[2])
window = 5000
ts_dir = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/1.sub_ts/"
)
out_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/delta_m/{model}_delta_m_bias_{bias}_sim.txt"

total_reps = defaultdict(dict)
for dir in ts_dir.glob("SEX~*/REC~*/MUT~5e-07"):
    REC = float(dir.parts[-2].split("~")[1])
    SEX = float(dir.parts[-3].split("~")[1])
    reps = []
    ts_files = sorted(dir.glob("*.trees"), key=lambda x: float(x.stem))
    for ts_file in tqdm.tqdm(ts_files):
        ts = tskit.load(str(ts_file))
        ts = ts.simplify(
            keep_input_roots=False
        )  # Otherwise: ValueError: Cannot rank trees with unary nodes
        if bias == 1:
            rep = calc_delta_m(ts)
        elif bias < 1:
            intervals = get_biased_intervals(ts, window, bias)
            ts = ts.keep_intervals(intervals)
            rep = calc_delta_m(ts)
        reps.append(rep)
    total_reps[SEX][REC] = reps

print(total_reps)

with open(out_file, "w") as f:
    json.dump(total_reps, f, sort_keys=True)
