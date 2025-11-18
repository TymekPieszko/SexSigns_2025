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

# python calc_tcomp_clustering_NEW.py MP 8 10 5000 0.1

model = sys.argv[1] 
n = int(sys.argv[2]) 
k = int(sys.argv[3])
window_size = int(sys.argv[4]) 
bias = float(sys.argv[5])
mut = 5e-07
N_anc = 1000
rec_anc = 1.0e-6
proc_num = n
rep_num = 100

target = Path(f"../sim_pipeline/sim_output/{model}/0.ts/")
out_file = Path(f"./tcomp/{model}_n_{n}_k_{k}_wind_{window_size}_bias_{bias}_clust.json")
checkpoint_file = Path(f"./tcomp/{model}_n_{n}_k_{k}_wind_{window_size}_bias_{bias}_clust.pkl")

if checkpoint_file.exists():
    with open(checkpoint_file, "rb") as f:
        total_reps = pickle.load(f)
else:
    total_reps = defaultdict(dict)

def calc_tcomp(ts, mut, k, n, bias, window_size, N_anc, rec_anc):
    seed = ph.seeds[int(ts.stem)]
    ts = tskit.load(str(ts))
    ts = pyslim.recapitate(ts, ancestral_Ne=N_anc, recombination_rate=rec_anc, random_seed=seed)
    tree_dict = {code: 0.0 for code in ph.codes}
    ts = ph.subsample_ts_seeded(ts, n, seed=seed)
    ts = msprime.sim_mutations(ts, rate=mut, random_seed=seed)
    for i, j in itertools.combinations(range(n), 2):
        ts_ij = ph.subsample_ts_to_given(ts, [i, j])
        intervals = ph.get_biased_intervals(ts_ij, window_size, bias)
        genos_ij, pos_ij = ph.get_biallelic_tcomp(ts_ij)
        for s, e in intervals:
            mask = (pos_ij > s) & (pos_ij < e)
            genos_ijk, pos_ijk = genos_ij[mask], pos_ij[mask]
            tree_dict = ph.run_clustering(genos_ijk, pos_ijk, tree_dict, k)
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
        func = partial(calc_tcomp, mut=mut, k=k, n=n, bias=bias, window_size=window_size, N_anc=N_anc, rec_anc=rec_anc)
        reps = list(tqdm(pool.imap(func, ts_files), total=len(ts_files)))
    total_reps[f"{sex}_{rec}"] = reps
    with open(checkpoint_file, "wb") as f:
        pickle.dump(total_reps, f)

print(total_reps)
with open(out_file, "w") as f:
    json.dump(total_reps, f)
