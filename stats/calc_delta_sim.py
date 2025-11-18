from pathlib import Path
import functions.phased as ph
import tskit, msprime, pyslim
from collections import defaultdict
import sys
import allel
import pickle, json
import itertools
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial
import warnings
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# python calc_delta_sim.py MP 8 5000 0.1 snps
# python calc_delta_sim.py MP 8 5000 1.0 snps
# python calc_delta_sim.py MP 8 5000 0.1 trees
# python calc_delta_sim.py MP 8 5000 1.0 trees

model = sys.argv[1] 
n = int(sys.argv[2]) 
window_size = int(sys.argv[3]) 
bias = float(sys.argv[4])
mode = sys.argv[5]
mut = 5e-07
N_anc = 1000
rec_anc = 1.0e-6
proc_num = n
rep_num = 100

target = Path(f"../sim_pipeline/sim_output/{model}/0.ts/")
out_file = Path(f"./delta/delta_{model}_n_{n}_wind_{window_size}_bias_{bias}_{mode}.json")
checkpoint_file = Path(f"./delta/delta_{model}_n_{n}_wind_{window_size}_bias_{bias}_{mode}.pkl")

if checkpoint_file.exists():
    with open(checkpoint_file, "rb") as f:
        total_reps = pickle.load(f)
else:
    total_reps = defaultdict(dict)

def calc_delta_master(ts, mut, n, bias, window_size, N_anc, rec_anc):
    seed = ph.seeds[int(ts.stem)]
    ts = tskit.load(str(ts))
    ts = pyslim.recapitate(ts, ancestral_Ne=N_anc, recombination_rate=rec_anc, random_seed=seed)
    ts = ph.subsample_ts_seeded(ts, n, seed=seed)
    ts = msprime.sim_mutations(ts, rate=mut, random_seed=seed)
    delta_pairs = []
    for i, j in itertools.combinations(range(n), 2):
        ts_ij = ph.subsample_ts_to_given(ts, [i, j])
        intervals = ph.get_biased_intervals(ts_ij, window_size, bias)
        # print(intervals)
        genos_ij, pos_ij = ph.get_biallelic_delta(ts_ij)
        delta_intervals = []
        for s, e in intervals:
            mask = (pos_ij > s) & (pos_ij < e)
            if sum(mask)==0:
                continue
            genos_ijk, pos_ijk = genos_ij[mask], pos_ij[mask]
            if mode == "snps":
                runs = ph.get_dbl_runs(genos_ijk, pos_ijk, pos_ijk[-1]+1)
            if mode == "trees":
                ts_ijk = ts_ij.keep_intervals([[s,e]])
                runs = ph.get_tree_runs(ts_ijk)
                # print(runs)
            delta = ph.run_delta(genos_ijk, pos_ijk, runs)
            delta_intervals.append(delta)
        # print(delta_intervals)
        delta_pairs.append(np.nanmean(delta_intervals))
    delta = np.nanmean(delta_pairs)   
    return delta

for dir in target.glob("SEX~*/REC~*"):
    sex = float(dir.parts[-2].split("~")[1])
    rec = float(dir.parts[-1].split("~")[1])
    if sex == 1.0:
        continue
    # if sex != 0.0:
    #     continue
    # if rec != 0.0:
    #     continue
    rec = float(dir.parts[-1].split("~")[1])
    if len(total_reps[f"{sex}_{rec}"]) == rep_num:
        continue
    ts_files = sorted(dir.glob("*.trees"), key=lambda x: float(x.stem))[:rep_num]
    with Pool(processes=proc_num) as pool:
        func = partial(calc_delta_master, mut=mut, n=n, bias=bias, window_size=window_size, N_anc=N_anc, rec_anc=rec_anc)
        reps = list(tqdm(pool.imap(func, ts_files), total=len(ts_files)))
    total_reps[f"{sex}_{rec}"] = reps
    with open(checkpoint_file, "wb") as f:
        pickle.dump(total_reps, f)

print(total_reps)
with open(out_file, "w") as f:
    json.dump(total_reps, f)
