import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# python plot_R2.py
file = "./out.txt"
df = pd.read_csv(file, sep="\t")
df = df.rename(index={0: "i"})
# df = df.sort_values(by="i", axis=1)
df = df.loc[:, ["Hi", "Fis", "R2", "CL", "AR", "SX", "CL_bias", "AR_bias", "SX_bias", "delta_m","dt_ratio", "old", "tcomp", "tcomp_bias", "treebased", "new"]] # The last 5 get reordered anyway.
# print(df)
df1 = df.iloc[:,:-5]
df2 = df.iloc[:,-5:].sort_values(by="i", axis=1)
df = pd.concat([df1,df2], axis=1)
print(df)
# # colors
pos = [0,1,2, 3.5,4.5,5.5, 7,8,9, 10.5,11.5, 13,14,15,16,17]
plt.figure(figsize=(8,5))
plt.bar(pos, df.iloc[0], width=0.8, edgecolor="black", color="white", hatch=[""]*3 + ["///"]*3 + ["/////"]*3 + ["'\\\\\\\\\\'"]*2 + ["/////", "", "'\\\\\\\\\\'", "xxxxx", "///"])
plt.xticks(ticks=pos, labels=df.columns, rotation=66)
plt.ylim(0,1)
plt.tight_layout()
plt.savefig("./R2.png", dpi=180)
plt.show()