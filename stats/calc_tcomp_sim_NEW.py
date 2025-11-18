from pathlib import Path
import sys
import tskit, pyslim, msprime
import pickle, json
from collections import defaultdict
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm
import functions.phased as ph
import itertools
from sexsigns_functions.calc import get_tree_dict_21, get_tree_composition_21, get_biased_intervals
import numpy as np
import warnings
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# python calc_tcomp_sim_NEW.py MP 8 5000 0.1
# python calc_tcomp_sim_NEW.py MP 8 5000 1.0
# THE BELOW IS FOR CHAPTER 3
# python calc_tcomp_sim_NEW.py MP 16 1000000 1.0
model = sys.argv[1]
n = int(sys.argv[2])
window_size = int(sys.argv[3])
bias = float(sys.argv[4])
mut = 5e-07
N_anc = 1000
rec_anc = 1.0e-6
rep_num = 100
# proc_num = n
proc_num = 32

target = Path(f"../sim_pipeline/sim_output/{model}/0.ts/")
out_file = Path(f"./tcomp/{model}_n_{n}_wind_{window_size}_bias_{bias}_sim.json")
checkpoint_file = Path(f"./tcomp/{model}_n_{n}_wind_{window_size}_bias_{bias}_sim.pkl")

if checkpoint_file.exists():
    with open(checkpoint_file, "rb") as f:
        total_reps = pickle.load(f)
else:
    total_reps = defaultdict(dict)   # total_reps[SEX][REC] = [tree_dicts]

def calc_tcomp_worker(ts_path, n, bias, window_size, N_anc, rec_anc):
    ts = tskit.load(str(ts_path))
    ts = pyslim.recapitate(ts, ancestral_Ne=N_anc, recombination_rate=rec_anc)
    # ts = ph.subsample_ts(ts, n)
    inds = np.random.choice(np.arange(N_anc), size=8, replace=True)
    tree_dict = {code: 0.0 for code in ph.codes}
    for i, j in itertools.combinations(inds, 2):
        ts_ij = ph.subsample_ts_to_given(ts, [i, j])
        ts_ij = ts_ij.simplify()
        if bias == 1.0:
            tree_dict = ph.calc_tcomp(ts_ij, tree_dict)
        if bias < 1.0:
            ts_ij = msprime.sim_mutations(ts_ij, rate=mut)
            intervals = ph.get_biased_intervals(ts_ij, window_size, bias)
            ts_ij = ts_ij.keep_intervals(intervals)
            tree_dict = ph.calc_tcomp(ts_ij, tree_dict)
    return tree_dict

for dir in target.glob("SEX~*/REC~*"):
    sex = float(dir.parts[-2].split("~")[1])
    if sex == 1.0:
        continue
    rec = float(dir.parts[-1].split("~")[1])
    if len(total_reps[f"{sex}_{rec}"]) == rep_num:
        continue
    ts_files = sorted(dir.glob("*.trees"), key=lambda x: float(x.stem))[:rep_num]
    with Pool(processes=proc_num) as pool:
        func = partial(calc_tcomp_worker, n=n, bias=bias, window_size=window_size, N_anc=N_anc, rec_anc=rec_anc)
        reps = list(tqdm(pool.imap(func, ts_files), total=len(ts_files)))
    total_reps[f"{sex}_{rec}"] = reps
    with open(checkpoint_file, "wb") as f:
        pickle.dump(total_reps, f)


with open(out_file, "w") as f:
    json.dump(total_reps, f)
