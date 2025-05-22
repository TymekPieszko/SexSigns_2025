from pathlib import Path
import tskit
import numpy as np
import pandas as pd
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import json

# python collect_stats.py MP
model = sys.argv[1]
stats_dir = Path(f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats")
out_file = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/collected/{model}_collected.tsv"
)
# 2inds - sample 2 inds per sim to caluclate Hi, Fis and LD_0-1000
# dv2 - use a different version of delta_m


in_files = [
    stats_dir / "Hi" / f"Hi_{model}_mut_5e-07_inds_100.txt",
    stats_dir / "Fis" / f"Fis_{model}_mut_5e-07_inds_100.txt",
    stats_dir / "LD_0-1000" / f"LD_0-1000_{model}_mut_5e-07_inds_100.txt",
    stats_dir / "delta_m" / f"{model}_delta_m_bias_1.0_sim.txt",
    stats_dir / "delta_m" / f"{model}_delta_m_bias_0.1_sim.txt",
    stats_dir / "delta_m" / f"{model}_delta_m_bias_0.1_SINGER.txt",
    stats_dir / "delta_m" / f"{model}_delta_m_bias_0.1_IQ_TREE.txt",
    stats_dir / "tree_composition" / f"{model}_sim_1.0.txt",
    stats_dir / "tree_composition" / f"{model}_sim_0.1.txt",
    stats_dir / "tree_composition" / f"{model}_SINGER_0.1.txt",
    stats_dir / "tree_composition" / f"{model}_IQ_TREE_0.1.txt",
]

# Parameter combinations
SEX_lst = [
    0.0,
    1.0e-05,
    3.162e-05,
    1.0e-04,
    3.162e-04,
    1.0e-03,
    3.162e-03,
    1.0e-02,
    3.162e-02,
    0.1,
    3.162e-01,
]
REC_lst = [0.0, 1.0e-09, 3.162e-09, 1.0e-08, 3.162e-08, 1.0e-07, 3.162e-07, 1.0e-06]
combos = [(SEX, REC) for SEX in SEX_lst for REC in REC_lst]


data_lst = []
for in_file in in_files:
    print(in_file)
    with open(in_file, "r") as f:
        data = json.load(f)
    f.close()
    if "tree_composition" in str(in_file):
        for tree in ["C", "A", "S"]:
            tmp = {
                SEX: {
                    REC: [calc_class_props(rep)[tree] for rep in reps]
                    for REC, reps in REC_dict.items()
                }
                for SEX, REC_dict in data.items()
            }
            data_lst.append(tmp)
    else:
        data_lst.append(data)
# print(len(data_lst))

df_total = pd.DataFrame()
for SEX, REC in combos:
    SEX = str(SEX)
    REC = str(REC)
    df_combo = pd.DataFrame(
        {
            "SEX": [SEX] * 100,
            "REC": [REC] * 100,
            "Hi": data_lst[0][SEX][REC],
            "Fis": data_lst[1][SEX][REC],
            "LD_0-1000": data_lst[2][SEX][REC],
            "D_1_0": data_lst[3][SEX][REC],
            "D_0_1": data_lst[4][SEX][REC],
            "D_SINGER": data_lst[5][SEX][REC],
            "D_IQ_TREE": data_lst[6][SEX][REC],
            "C_1_0": data_lst[7][SEX][REC],
            "A_1_0": data_lst[8][SEX][REC],
            "S_1_0": data_lst[9][SEX][REC],
            "C_0_1": data_lst[10][SEX][REC],
            "A_0_1": data_lst[11][SEX][REC],
            "S_0_1": data_lst[12][SEX][REC],
            "C_SINGER": data_lst[13][SEX][REC],
            "A_SINGER": data_lst[14][SEX][REC],
            "S_SINGER": data_lst[15][SEX][REC],
            "C_IQ_TREE": data_lst[16][SEX][REC],
            "A_IQ_TREE": data_lst[17][SEX][REC],
            "S_IQ_TREE": data_lst[18][SEX][REC],
        }
    )
    print(df_combo)
    df_total = pd.concat([df_total, df_combo], axis=0)

df_total.to_csv(out_file, sep="\t", index=False)
