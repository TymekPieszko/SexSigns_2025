from pathlib import Path
from functions import *
import time

fasta_dir = Path(snakemake.input[0])
tree_dir = Path(snakemake.output[0])
threads = snakemake.threads

start_time = time.time()
# i = 0
tree_dir.mkdir(exist_ok=True)
for fasta_f in fasta_dir.glob("*.fa.gz"):
    # if i > 9:
    #     break
    infer_tree(fasta_f, tree_dir, threads)
    for file in tree_dir.iterdir():
        if file.suffix != ".treefile":
            file.unlink()
#     i += 1
end_time = time.time()
print(f"Elapsed time: {end_time - start_time}")
