from pathlib import Path
import tskit
import numpy as np
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import json

# python get_tree_composition_sim.py MP 1.0
# python get_tree_composition_sim.py CF 1.0
model = sys.argv[1]
bias = float(sys.argv[2])
window = 5000

ts_dir = Path(
    f"../sim_pipeline/sim_output/{model}/1.sub_ts/"
)
out_file = f"./tree_composition/{model}_sim_{bias}.txt"


total_dict = defaultdict(dict)
for dir in ts_dir.glob("SEX~*/REC~*/MUT~*"):
    SEX = float(dir.parts[-3].split("~")[1])
    REC = float(dir.parts[-2].split("~")[1])
    # print(SEX, REC)
    reps = []
    ts_files = sorted(dir.glob("*.trees"), key=lambda x: float(x.stem))
    for ts_file in tqdm.tqdm(ts_files):
        tree_dict = get_tree_dict_21()
        ts = tskit.load(str(ts_file))
        if bias == 1:
            tree_dict = get_tree_composition_21(ts, tree_dict)
        elif bias < 1:
            intervals = get_biased_intervals(ts, window, bias)
            ts = ts.keep_intervals(intervals)
            tree_dict = get_tree_composition_21(ts, tree_dict)
        reps.append(tree_dict)
    total_dict[SEX][REC] = reps

with open(out_file, "w") as f:
    json.dump(total_dict, f)
