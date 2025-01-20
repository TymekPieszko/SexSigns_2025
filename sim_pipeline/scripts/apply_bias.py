import numpy as np
import tskit
import os, sys
from functions_smk import *

ts = snakemake.input[0]
ts_dir = snakemake.output[0]  # Output directory to store phased VCFs
BIAS = float(snakemake.wildcards.BIAS)
WINDOW_LEN = int(snakemake.params.WINDOW_LEN)

ts = tskit.load(str(ts))
windows = np.arange(0, ts.sequence_length + WINDOW_LEN, WINDOW_LEN)
sample_sets = [[0, 1], [2, 3]]


het_windows = ts.diversity(
    sample_sets=sample_sets, windows=windows, mode="site", span_normalise=True
)
het_windows = np.mean(
    het_windows, axis=1
)  # Calc mean H_I per individual; take top BIAS proportion of windows

index = int(len(het_windows) * BIAS)
het_filter = np.argsort(het_windows)[-index:]
intervals = [[windows[i], windows[i + 1]] for i in range(len(windows) - 1)]
intervals = np.array(intervals)[het_filter]

os.makedirs(ts_dir, exist_ok=True)
for int, i in zip(intervals, het_filter):
    int_ts = ts.keep_intervals([int]).trim()
    int_ts.dump(os.path.join(ts_dir, f"{i}.trees"))
