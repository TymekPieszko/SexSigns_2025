from pathlib import Path
import newick
import tskit
import numpy as np
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import pickle, json

# python calc_delta_m_IQ_TREE.py MP 0.1
model = sys.argv[1]
bias = float(sys.argv[2])

iqtree_dir = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/3.2.IQ_TREE/"
)
fasta_dir = f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/2.2.fasta/"
checkpoint_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/delta_m/{model}_delta_m_bias_{bias}_IQ_TREE.pkl"
out_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/delta_m/{model}_delta_m_bias_{bias}_IQ_TREE.txt"


if Path(checkpoint_file).exists():
    with open(checkpoint_file, "rb") as f:
        total_dict = pickle.load(f)
else:
    total_dict = defaultdict(dict)
for param_dir in iqtree_dir.glob(f"SEX~*/REC~*/MUT~*"):
    SEX = float(param_dir.parts[-3].split("~")[1])
    REC = float(param_dir.parts[-2].split("~")[1])
    if (
        SEX in total_dict
        and REC in total_dict[SEX]
        and len(total_dict[SEX][REC]) == 100
    ):
        print(f"Skipping {SEX} + {REC}")
        continue
    reps = []
    rep_dirs = sorted(param_dir.glob("*/sub"), key=lambda x: float(x.parts[-2]))
    for rep_dir in tqdm.tqdm(rep_dirs):
        # print(rep_dir)
        i = int(rep_dir.parts[-2])
        # if i > 1:
        #     continue
        dm_lst = []
        hi_lst = []
        for newick_f in rep_dir.glob("*.treefile"):
            # print(newick_f)
            # if int(newick_f.stem.split("_")[1]) % 100000 != 0:
            #     continue

            ###############################
            ### Calculate delta_m
            tree = newick.loads(newick_f.read_text())[0]
            tree = rename_nodes(tree)
            # print(tree.ascii_art())
            # print(newick.dumps(tree))
            dm = calc_delta_m_iq(tree)
            dm_lst.append(dm)

            ###############################
            ### Calculate Hi
            fasta_f = (
                fasta_dir + f"SEX~{SEX}/REC~{REC}/MUT~5e-07/{i}/{newick_f.stem}.fa.gz"
            )
            samples = newick.loads(newick_f.read_text())[
                0
            ].get_leaf_names()  # Need to load a fresh newick
            # print(samples)
            hi = calc_hi_from_fasta(fasta_f, samples)
            hi_lst.append(hi)
        # print(np.sum(np.isnan(dm_lst)), np.sum(np.isnan(hi_lst)))
        dm_lst = filter_by_hi(dm_lst, hi_lst, bias)
        # print(dm_lst)
        reps.append(np.nanmean(dm_lst))
    total_dict[SEX][REC] = reps

    with open(checkpoint_file, "wb") as f:
        pickle.dump(total_dict, f)

with open(out_file, "w") as f:
    json.dump(total_dict, f, sort_keys=True)
