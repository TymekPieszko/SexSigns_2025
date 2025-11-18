import sys, json
import pandas as pd
import functions.phased as ph
from functions.plot import sex_lst, rec_lst, combo_keys
from sexsigns_functions.plot import calc_class_props

def flatten_dict(d):
    out = {}
    for k1, sub in d.items():
        for k2, v in sub.items():
            out[f"{k1}_{k2}"] = v
    return out

# python collect_stats.py MP
model = str(sys.argv[1])
# tcompTotal_file = f"../stats/tcomp/{model}_n_8_k_10_wind_5000_bias_1.0_clust.json"
# tcompBias_file = f"../stats/tcomp/{model}_n_8_k_10_wind_5000_bias_0.1_clust.json"
tcompTotal_file = f"../stats/tcomp/{model}_n_8_wind_5000_bias_1.0_sim.json"
tcompBias_file = f"../stats/tcomp/{model}_IQ_TREE_0.1.json"
delta_m_file = f"../stats/delta/delta_{model}_n_8_wind_5000_bias_0.1_snps.json"
dt_ratio_file = f"../stats/kappa/kappa_{model}_inds_8.json"
hi_file = f"../stats/Hi/Hi_{model}_mut_5e-07_inds_8.json"
fis_file = f"../stats/Fis/Fis_{model}_mut_5e-07_inds_8.json"
r2_file = f"../stats/LD_0-1000/LD_0-1000_{model}_mut_5e-07_inds_8.json"
out_file = f"./{model}_stats.tsv"

trees = ["CL", "AR", "SX"]
with open(tcompTotal_file, "r") as f:
    tcomp_data = json.load(f)
CL_data_1 = {k: [ph.calc_category_props(tree_dict)[trees[0]] for tree_dict in tcomp_data[k]] for k in combo_keys}
AR_data_1 = {k: [ph.calc_category_props(tree_dict)[trees[1]] for tree_dict in tcomp_data[k]] for k in combo_keys}
SX_data_1 = {k: [ph.calc_category_props(tree_dict)[trees[2]] for tree_dict in tcomp_data[k]] for k in combo_keys}
# tcomp_data = flatten_dict(tcomp_data)
# CL_data_1 = {k: [calc_class_props(tree_dict)[trees[0]] for tree_dict in tcomp_data[k]] for k in combo_keys}
# AR_data_1 = {k: [calc_class_props(tree_dict)[trees[1]] for tree_dict in tcomp_data[k]] for k in combo_keys}
# SX_data_1 = {k: [calc_class_props(tree_dict)[trees[2]] for tree_dict in tcomp_data[k]] for k in combo_keys}
with open(tcompBias_file, "r") as f:
    tcomp_data = json.load(f)
trees = ["C", "A", "S"]
# CL_data_2 = {k: [ph.calc_category_props(tree_dict)[trees[0]] for tree_dict in tcomp_data[k]] for k in combo_keys}
# AR_data_2 = {k: [ph.calc_category_props(tree_dict)[trees[1]] for tree_dict in tcomp_data[k]] for k in combo_keys}
# SX_data_2 = {k: [ph.calc_category_props(tree_dict)[trees[2]] for tree_dict in tcomp_data[k]] for k in combo_keys}
tcomp_data = flatten_dict(tcomp_data)
CL_data_2 = {k: [calc_class_props(tree_dict)[trees[0]] for tree_dict in tcomp_data[k]] for k in combo_keys}
AR_data_2 = {k: [calc_class_props(tree_dict)[trees[1]] for tree_dict in tcomp_data[k]] for k in combo_keys}
SX_data_2 = {k: [calc_class_props(tree_dict)[trees[2]] for tree_dict in tcomp_data[k]] for k in combo_keys}

with open(delta_m_file, "r") as f:
    delta_m_data = json.load(f)
with open(dt_ratio_file, "r") as f:
    dt_ratio_data = json.load(f)
with open(hi_file, "r") as f:
    hi_data = json.load(f)
with open(fis_file, "r") as f:
    fis_data = json.load(f)
with open(r2_file, "r") as f:
    r2_data = json.load(f)

table = []
for sex in sex_lst:
    for rec in rec_lst:
        a = CL_data_1[f"{sex}_{rec}"]
        b = AR_data_1[f"{sex}_{rec}"]
        c = SX_data_1[f"{sex}_{rec}"]
        d = CL_data_2[f"{sex}_{rec}"]
        e = AR_data_2[f"{sex}_{rec}"]
        f = SX_data_2[f"{sex}_{rec}"]
        g = delta_m_data[f"{sex}_{rec}"]
        h = dt_ratio_data[f"{sex}_{rec}"]
        i = hi_data[f"{sex}_{rec}"]
        j = fis_data[f"{sex}_{rec}"]
        k = r2_data[f"{sex}_{rec}"]
        for idx in range(100):
            table.append([sex, rec, a[idx], b[idx], c[idx], d[idx], e[idx], f[idx], g[idx], h[idx], i[idx], j[idx], k[idx]])

cols = ['sex', 'rec', 'CL', 'AR', 'SX', 'CL_bias', 'AR_bias', 'SX_bias', 'delta_m', 'dt_ratio', 'Hi', 'Fis', 'R2']
df = pd.DataFrame(table, columns=cols)
df.to_csv(out_file, sep='\t', index=False)
