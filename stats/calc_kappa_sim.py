from pathlib import Path
import tskit, msprime, pyslim
import numpy as np
import allel
import tqdm
import functions.unphased as uph
import sys
from collections import defaultdict
import json
from multiprocessing import Pool
from functools import partial
import warnings
warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)

# python calc_kappa_sim.py MP 8
# python calc_kappa_sim.py CF 8
model = sys.argv[1]
n = int(sys.argv[2])
mut = 5e-07
N_anc = 1000
rec_anc = 1.0e-6
proc_num = n  # number of worker processes

ts_dir = Path(
    f"../sim_pipeline/sim_output/{model}/0.ts/"
)
out_file = f"./kappa/kappa_{model}_inds_{n}.json"

total_reps = defaultdict(dict)

def calc_dt_ratio(ts_file, sex, rec, mut, n, N_anc, rec_anc):
    seed = uph.seeds[int(ts_file.stem)]
    ts = tskit.load(str(ts_file))
    ts = pyslim.recapitate(ts, ancestral_Ne=N_anc, recombination_rate=rec_anc, random_seed=seed)
    ts = uph.subsample_ts_seeded(ts, n, seed=seed)
    ts = msprime.sim_mutations(ts, rate=mut, random_seed=seed)
    genos = allel.HaplotypeArray(ts.genotype_matrix()).to_genotypes(ploidy=2)
    arr = uph.calc_dt_ratio(genos, n)
    rep = np.nanmean(arr[~np.eye(n, dtype=bool)])
    return rep

for dir in ts_dir.glob("SEX~*/REC~*"):
    sex = float(dir.parts[-2].split("~")[1])
    if sex == 1.0:
        continue
    rec = float(dir.parts[-1].split("~")[1])
    print(sex, rec)

    ts_files = sorted(dir.glob("*.trees"), key=lambda x: float(x.stem))

    with Pool(processes=proc_num) as pool:
        func = partial(calc_dt_ratio, sex=sex, rec=rec, mut=mut,
                       n=n, N_anc=N_anc, rec_anc=rec_anc)
        reps = list(tqdm.tqdm(pool.imap(func, ts_files), total=len(ts_files)))

    total_reps[f"{sex}_{rec}"] = reps

with open(out_file, "w") as f:
    json.dump(total_reps, f, sort_keys=True)
