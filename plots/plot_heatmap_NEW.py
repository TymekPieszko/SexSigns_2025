from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime
import json, sys
from functions.plot import plot_heatmap, params, combo_keys, sex_lst, rec_lst


# python plot_heatmap_NEW.py kappa kappa_CF_inds_8
# python plot_heatmap_NEW.py delta delta_MP_n_8_wind_5000_bias_0.1_snps
# python plot_heatmap_NEW.py delta delta_MP_bias_0.1_SINGER

stat = sys.argv[1]
target = sys.argv[2]
in_file = f"../stats/{stat}/{target}.json"
plot_file = Path(f"./heatmaps/{target}.png")

with open(in_file, "r") as f:
    data = json.load(f)
print(data)
data = {k: np.nanmean(data[k]) for k in combo_keys}
data = {
    float(r): {
        float(s): data.get(f"{s}_{r}", np.nan)
        for s in sex_lst
    }
    for r in rec_lst
}
df = pd.DataFrame(data).T.sort_index(ascending=False)
print(df)

title = f"{target}; {datetime.today().strftime('%Y-%m-%d %H:%M:%S')}"
plot_heatmap(
    df=df,
    figsize=(14, 10),
    cmap="plasma",
    plot_file=plot_file,
    title=title,
    xlab=r"$\mathit{\sigma}$",
    ylab=r"$\mathit{\gamma}$",
    label_font=66,
    tick_font=40,
    title_font=20,
    vmin=np.min(df),
    vmax=np.max(df),
    annot=False,
)
