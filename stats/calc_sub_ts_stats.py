from pathlib import Path
from collections import defaultdict
import tskit
import pandas as pd
import tqdm
import numpy as np

model = "MP"
out_dir = "./sub_ts_stats/"
total_blocks = defaultdict(dict)
total_SNVs = defaultdict(dict)
for dir in Path(f"../sim_pipeline/sim_output/{model}/1.sub_ts").glob(
    "SEX~*/REC~*/MUT~*/"
):
    SEX = float(dir.parts[-3].split("~")[1])
    REC = float(dir.parts[-2].split("~")[1])
    print(SEX, REC)
    blocks_ts = []
    SNVs_ts = []
    for ts in tqdm.tqdm(dir.glob("*trees")):
        # if float(ts.stem) > 10:
        #     continue
        ts = tskit.load(ts)
        blocks_ts.append(ts.num_trees)
        SNVs_ts.append(ts.num_sites)
    print(len(blocks_ts), len(SNVs_ts))
    total_blocks[SEX][REC] = np.nanmean(blocks_ts)
    total_SNVs[SEX][REC] = np.nanmean(SNVs_ts)
# print(total_blocks)
# print(total_SNVs)
total_blocks = (
    pd.DataFrame(total_blocks).sort_index().sort_index(axis=1).iloc[:, :-1][::-1]
)
total_SNVs = pd.DataFrame(total_SNVs).sort_index().sort_index(axis=1).iloc[:, :-1][::-1]
total_blocks.to_csv(out_dir + f"{model}_blocks.txt", sep="\t", header=True, index=True)
total_SNVs.to_csv(out_dir + f"{model}_SNVs.txt", sep="\t", header=True, index=True)
