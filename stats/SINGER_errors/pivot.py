import pandas as pd

model = "CF"
target = f"./{model}_errors.txt"
df = pd.read_csv(target, sep="\t", header=None)
df = df.pivot(index=1, columns=0, values=2).iloc[:, :-1][::-1]
df.to_csv(target.replace(".txt", "_pivot.txt"), sep="\t", header=True, index=True)
