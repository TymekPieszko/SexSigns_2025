import tskit
from pathlib import Path
from functions import *

ts = Path(snakemake.input.ts)
vcf = Path(snakemake.input.vcf)
ts_dir = Path(snakemake.output.ts_dir)
arg_dir = Path(str(ts_dir).replace("ts", "arg"))

arg_log = Path(snakemake.log.arg_log)
ts_log = Path(snakemake.log.ts_log)

L = int(snakemake.params.L)
WINDOW_LEN = int(snakemake.params.WINDOW_LEN)
BIAS = float(snakemake.wildcards.BIAS)
ts = tskit.load(str(ts))
intervals = get_biased_intervals(ts, WINDOW_LEN, BIAS)
# print(intervals)

# Infer ARG params
singer_path = "singer_master"
Ne = "1000" # 3.1.SINGER_OLD: Ne = "100000"
m = "5.0e-7" # 3.1.SINGER_OLD: m = "1.0e-8"
n = "100"
thin = "20"

# Convert to ts params
convert_to_tskit_path = "convert_to_tskit"  # "./convert_to_tskit.py"
start_index = "50"
end_index = "100"
step = "5"

for interval in intervals:
    vcf_prefix = str(vcf.parent / vcf.stem)
    name = f"{int(interval[0])}_{int(interval[1])}"
    arg_dir_i = arg_dir / name
    arg_prefix = str(arg_dir_i / name)
    arg_dir_i.mkdir(parents=True, exist_ok=True)

    ts_dir_i = ts_dir / name
    ts_prefix = str(ts_dir_i / name)
    ts_dir_i.mkdir(parents=True, exist_ok=True)

    start_pos = str(interval[0])
    end_pos = str(interval[1])

    with open(arg_log, "a") as log:
        run_singer(
            singer_path=singer_path,
            Ne=Ne,
            m=m,
            vcf_prefix=vcf_prefix,
            arg_prefix=arg_prefix,
            start=start_pos,
            end=end_pos,
            n=n,
            thin=thin,
            log=log,
        )
    with open(ts_log, "a") as log:
        run_convert_to_tskit(
            convert_to_tskit_path=convert_to_tskit_path,
            arg_prefix=arg_prefix,
            ts_prefix=ts_prefix,
            start=start_index,
            end=end_index,
            step=step,
            log=log,
        )
