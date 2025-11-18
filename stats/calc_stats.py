from pathlib import Path
import tskit, msprime, pyslim
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import json
import allel
import numpy as np
from multiprocessing import Pool
import warnings
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# python calc_stats.py MP Hi 100
# python calc_stats.py MP Fis 100
# python calc_stats.py MP LD_0-1000 100
model = sys.argv[1]
stat = sys.argv[2]
NUM_INDS = int(sys.argv[3])
N_ANC = 1000
REC_ANC = 1.0e-6
MUT = 5e-07
ts_dir = Path(f"../sim_pipeline/sim_output/{model}/0.ts/")
out_file = f"./{stat}/{stat}_{model}_mut_{MUT}_inds_{NUM_INDS}.json"


def calc_r2_allel(ts, K, d_max):
    genos = allel.HaplotypeArray(ts.genotype_matrix()).to_genotypes(ploidy=2)
    pos = ts.sites_position
    n_var = len(genos)
    K = min(K, n_var)  # cannot sample more targets than exist
    r2_total = []
    for _ in range(K):
        target = np.random.randint(0, n_var)
        target_geno = genos[target]
        target_pos = pos[target]
        dist = np.abs(pos - target_pos)
        valid = (dist > 0) & (dist < d_max)
        if not np.any(valid):
            continue

        # make 2D arrays
        target_geno_n = target_geno.to_n_alt()[None, :]    # shape (1, n_samples)
        other_genos_n = genos[valid].to_n_alt()            # shape (m, n_samples)

        r = allel.rogers_huff_r_between(target_geno_n, other_genos_n)
        r2_total.extend(np.ravel(r) ** 2)

    if len(r2_total) == 0:
        return np.nan
    else:
        return float(np.nanmean(r2_total))


def process_ts(args):
    """Worker: load ts_file, compute the requested stat."""
    ts_file, stat = args
    ts = tskit.load(str(ts_file))
    ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC_ANC)
    ts = subsample_ts(ts, NUM_INDS)
    ts = msprime.sim_mutations(ts, rate=MUT)

    if stat == "Hi":
        return calc_hi(ts)
    elif stat == "Fis":
        return calc_fis(ts)
    elif stat == "LD_0-1000":
        return calc_r2_allel(ts, 500, 1000)
    # elif stat == "LD_9000-10000":
    #     return calc_r2(ts, 9000, 10000)
    else:
        return None


total_reps = defaultdict(dict)

for dir in ts_dir.glob("SEX~*/REC~*"):
    rec = float(dir.parts[-1].split("~")[1])
    sex = float(dir.parts[-2].split("~")[1])
    ts_files = sorted(dir.glob("*.trees"), key=lambda x: float(x.stem))

    with Pool() as pool:
        reps = list(tqdm.tqdm(pool.imap(process_ts, [(f, stat) for f in ts_files]),
                              total=len(ts_files)))
    total_reps[f"{sex}_{rec}"] = reps

print(total_reps)
with open(out_file, "w") as f:
    json.dump(total_reps, f, sort_keys=True)
