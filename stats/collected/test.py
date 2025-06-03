import pandas as pd
import numpy as np

data = "./MP_collected.tsv"
df = pd.read_csv(data, sep="\t", header=0)
print(df.head())
small = df.loc[(df["SEX"] == 3.162e-5) & (df["REC"] == 0.0)]
print(len(small.loc[small["S_1_0"] == 0.0]))
