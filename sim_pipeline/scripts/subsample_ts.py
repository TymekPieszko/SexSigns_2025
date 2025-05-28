import tskit, pyslim, msprime, os, gzip
from pathlib import Path
import numpy as np
from functions import *

### Note
# Ideally, you should probably (1) overlay the full ts, (2) get 100-ind ts, then (3) get 2-ind ts from 100-ind. But order does not matter.

### Get paths / params
ts_in = Path(snakemake.input[0])
ts_out = snakemake.output[0]
vcf_out = snakemake.output[1]
fasta_dir_out = snakemake.output[2]
MUT = float(snakemake.wildcards.MUT)
N_ANC = int(snakemake.params.N_ANC)
L = float(snakemake.params.L)
REC = 1 / L
WINDOW_LEN = int(snakemake.params.WINDOW_LEN)
windows = np.arange(0, L + WINDOW_LEN, WINDOW_LEN)

### Subsample the ts, output a sub-ts + vcf
ts = tskit.load(ts_in)
ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC)
ts2, inds2 = subsample_ts(ts, 2)
ts2 = msprime.sim_mutations(ts2, rate=MUT)
# Output a 2-ind ts
ts2.dump(ts_out)
# Output a vcf for the 2-ind ts
with open(vcf_out, "w") as f:
    ts2.write_vcf(f, allow_position_zero=True)
f.close()


### Output window-wise fasta containing the 2 target inds
inds100 = np.concatenate(
    [
        inds2,
        np.random.choice(np.setdiff1d(np.arange(1000), inds2), 98, replace=False),
    ]
)
indices = np.argsort(np.argsort(inds100))[0:2]
# Record the indices of target inds for extracting them from the 100-ind tree inferred with IQ-TREE
os.makedirs(fasta_dir_out, exist_ok=True)  # Same dir as fastas; see below
### BUT I THINK THIS STEP IS NOT NECESSARY!!!!!!!!!! Snakemake creates directories for you!!!
with open(fasta_dir_out + "/inds.txt", "w") as f:
    for i in indices:
        f.write(f"{i}\t")
f.close()

ts100 = subsample_ts_to_given(ts, inds100)
ts100 = msprime.sim_mutations(ts100, rate=MUT)
for i, j in zip(windows[:-1], windows[1:]):
    ts100_window = ts100.keep_intervals([[i, j]]).trim()
    with gzip.open(fasta_dir_out + f"/{int(i)}_{int(j)}.fa.gz", "wt") as f:
        ts100_window.write_fasta(f)


###################################
### EMERGENCY VERSION 2025-04-29
# uncomment all lines at once
###################################

# import tskit, pyslim, msprime, os, gzip
# from pathlib import Path
# import numpy as np
# from functions import *

# ### Note
# # Ideally, you should probably (1) overlay the full ts, (2) get 100-ind ts, (3) get 2-ind ts from 100-ind, etc.. However, you don't want to rerun SINGER at this point (2025-04-29). Does not seem to be a big deal as well.

# ### Get paths / params
# ts_in = Path(snakemake.input[0])
# # ts_out = snakemake.output[0]
# # vcf_out = snakemake.output[1]
# fasta_dir_out = snakemake.output[0]
# MUT = float(snakemake.wildcards.MUT)
# N_ANC = int(snakemake.params.N_ANC)
# L = float(snakemake.params.L)
# REC = 1 / L
# WINDOW_LEN = int(snakemake.params.WINDOW_LEN)
# windows = np.arange(0, L + WINDOW_LEN, WINDOW_LEN)

# #######################################
# ### FOR NOW, DON'T RUN
# ### Load and subsample the ts
# # ts = tskit.load(ts_in)
# # ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC)
# # ts2, inds2 = subsample_ts(ts, 2)
# # ts2 = msprime.sim_mutations(ts2, rate=MUT)
# # # Output a 2-ind ts
# # ts2.dump(ts_out)
# # # Output a vcf for the 2-ind ts
# # with open(vcf_out, "w") as f:
# #     ts2.write_vcf(f, allow_position_zero=True)
# # Write the inds to file
# os.makedirs(fasta_dir_out, exist_ok=True)  # Same dir as fastas; see below
# # with open(fasta_dir_out + "/inds.txt", "w") as f:
# #     for i in inds2:
# #         f.write(f"{i}\t")
# #######################################

# SEX = float(ts_in.parts[-3].split("~")[1])
# REC = float(ts_in.parts[-2].split("~")[1])
# print(SEX, REC)
# sub_ts = f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/CF/1.sub_ts/SEX~{SEX}/REC~{REC}/MUT~5e-07/{ts_in.name}"
# print(sub_ts)
# ts = tskit.load(ts_in)
# ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC)
# ts2 = tskit.load(sub_ts)
# inds2 = []
# ts_inds = list(ts.individuals())

# for i in range(2):
#     meta_i = ts2.individual(i).metadata
#     for j in range(1000):
#         meta_j = ts_inds[j].metadata
#         if meta_i == meta_j:
#             inds2.append(j)
#         # print(inds2)

# inds100 = np.concatenate(
#     [
#         inds2,
#         np.random.choice(np.setdiff1d(np.arange(1000), inds2), 98, replace=False),
#     ]
# )
# inds100 = [int(i) for i in inds100]
# indices = np.argsort(np.argsort(inds100))[0:2]
# with open(fasta_dir_out + "/inds.txt", "w") as f:
#     for i in indices:
#         f.write(f"{i}\t")

# f.close()

# ### Output per-window fastas for a 100-ind ts
# # Include the 2 previously sampled inds!!!
# # inds100 = np.concatenate(
# #     [inds2, np.random.choice(np.setdiff1d(np.arange(1000), inds2), 98, replace=False)]
# # )
# print(inds100)
# ts100 = subsample_ts_to_given(ts, inds100)
# ts100 = msprime.sim_mutations(ts100, rate=MUT)
# for i, j in zip(windows[:-1], windows[1:]):
#     ts100_window = ts100.keep_intervals([[i, j]]).trim()
#     with gzip.open(fasta_dir_out + f"/{int(i)}_{int(j)}.fa.gz", "wt") as f:
#         ts100_window.write_fasta(f)
