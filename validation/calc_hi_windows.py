from pathlib import Path
import tskit, msprime, pyslim
import tqdm
from sexsigns_functions.calc import subsample_ts, calc_hi_windows
import sys
from collections import defaultdict
import json

# python calc_hi_windows.py MP 100
# python calc_hi_windows.py CF 100
model = sys.argv[1]
NUM_INDS = int(sys.argv[2])
N_ANC = 1000
REC_ANC = 1.0e-6
MUT = 5e-07
window_length = 10000


ts_dir = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/validation/validation_sims/sim_output/{model}/0.ts/"
)
# ts_dir = Path(
#     f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/0.ts/"
# )
out_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/validation/Hi_windows_{model}_mut_{MUT}_inds_{NUM_INDS}.txt"

total_reps = defaultdict(dict)
for dir in ts_dir.glob("SEX~*/REC~*"):
    REC = float(dir.parts[-1].split("~")[1])
    SEX = float(dir.parts[-2].split("~")[1])
    # print(SEX, REC)
    reps = []
    for ts_file in tqdm.tqdm(dir.glob("*.trees")):
        # if float(ts_file.stem) > 30:
        #     continue
        ts = tskit.load(str(ts_file))
        ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC_ANC)
        ts = subsample_ts(ts, NUM_INDS)
        ts = ts.simplify()
        ts = msprime.sim_mutations(ts, rate=MUT)
        hi_windows = calc_hi_windows(ts, window_length)
        reps.append(
            list(hi_windows)
        )  # TypeError: Object of type ndarray is not JSON serializable
    total_reps[SEX][REC] = reps

with open(out_file, "w") as f:
    json.dump(total_reps, f, sort_keys=True)
