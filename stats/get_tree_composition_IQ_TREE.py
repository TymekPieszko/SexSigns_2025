from pathlib import Path
import newick
import numpy as np
import tqdm
from sexsigns_functions.calc import *
import sys
from collections import defaultdict
import pickle, json

# python get_tree_composition_IQ_TREE.py MP 0.1
model = sys.argv[1]
bias = float(sys.argv[2])
window = 5000

iqtree_dir = Path(
    f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/3.2.IQ_TREE/"
)
ts_dir = f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/1.sub_ts/"
fasta_dir = f"/data/biol-bdelloids/scro4331/SexSigns_2025/sim_pipeline/sim_output/{model}/2.2.fasta/"
checkpoint_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/tree_composition/{model}_IQ_TREE_{bias}.pkl"
out_file = f"/data/biol-bdelloids/scro4331/SexSigns_2025/stats/tree_composition/{model}_IQ_TREE_{bias}.txt"


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
        cat_dicts = []
        hi_values = []
        for newick_f in rep_dir.glob("*.treefile"):
            # print(newick_f)
            # if int(newick_f.stem.split("_")[1]) % 100000 != 0:
            #     continue
            # CAT DICT
            tree = newick.loads(newick_f.read_text())[0]
            # print(tree.ascii_art())
            cat_dict = get_tree_category(tree, window)
            # print(cat_dict)
            cat_dicts.append(cat_dict)

            # HI
            fasta_f = (
                fasta_dir + f"SEX~{SEX}/REC~{REC}/MUT~5e-07/{i}/{newick_f.stem}.fa.gz"
            )
            samples = newick.loads(newick_f.read_text())[
                0
            ].get_leaf_names()  # Need to load a fresh newick
            hi = calc_hi_from_fasta(fasta_f, samples)
            hi_values.append(hi)
        cat_dicts = filter_by_hi(cat_dicts, hi_values, bias)
        tree_dict = sum_cat_dicts(cat_dicts)
        reps.append(tree_dict)
    total_dict[SEX][REC] = reps

    with open(checkpoint_file, "wb") as f:
        pickle.dump(total_dict, f)

with open(out_file, "w") as f:
    json.dump(total_dict, f, sort_keys=True)
