from pathlib import Path
import tskit, msprime, pyslim
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import json

# python calc_Fis_analysis.py MP
model = sys.argv[1]
N_ANC = 1000
REC_ANC = 1.0e-6
INDS_lst = [5, 10, 15]
MUT_lst = [5e-09, 5e-08, 5e-07]
combos = [(i, j) for i in INDS_lst for j in MUT_lst]
ts_dir = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/0.ts/"
)

for combo in combos:
    INDS = combo[0]
    MUT = combo[1]
    out_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/Fis/Fis_{model}_mut_{MUT}_inds_{INDS}.txt"
    total_reps = defaultdict(dict)
    for dir in ts_dir.glob("SEX~*/REC~*"):
        REC = float(dir.parts[-1].split("~")[1])
        SEX = float(dir.parts[-2].split("~")[1])
        if SEX != 0.0:
            continue
        reps = []
        for ts_file in tqdm.tqdm(dir.glob("*.trees")):
            ts = tskit.load(str(ts_file))
            ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC_ANC)
            ts = subsample_ts(ts, INDS)
            ts = msprime.sim_mutations(ts, rate=MUT)
            rep = calc_fis(ts)
            reps.append(rep)
        total_reps[SEX][REC] = reps
    with open(out_file, "w") as f:
        json.dump(total_reps, f, sort_keys=True)
