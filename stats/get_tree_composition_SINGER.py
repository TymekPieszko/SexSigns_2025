from pathlib import Path
import tskit
import numpy as np
import tqdm
from calc_functions import *
import sys
from collections import defaultdict
import json

# python get_tree_composition_SINGER.py MP 0.1
model = sys.argv[1]
bias = float(sys.argv[2])
window = 5000

ts_dir = Path(
    f"../sim_pipeline/sim_output/{model}/3.1.SINGER/"
)
out_file = f"./tree_composition/{model}_SINGER_{bias}.txt"


total_dict = defaultdict(dict)
for param_dir in ts_dir.glob(f"SEX~*/REC~*/MUT~*/BIAS~{bias}"):
    SEX = float(param_dir.parts[-4].split("~")[1])
    REC = float(param_dir.parts[-3].split("~")[1])
    # print(SEX, REC)
    reps = []
    rep_dirs = sorted(param_dir.glob("*/ts"), key=lambda x: float(x.parts[-2]))
    for rep_dir in tqdm.tqdm(rep_dirs):
        tree_dict = get_tree_dict_21()
        for ts in rep_dir.rglob("*.trees"):
            ts = tskit.load(str(ts))
            tree_dict = get_tree_composition_21(ts, tree_dict)
        reps.append(tree_dict)
    total_dict[SEX][REC] = reps

with open(out_file, "w") as f:
    json.dump(total_dict, f)
