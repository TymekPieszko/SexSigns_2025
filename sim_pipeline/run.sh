# module load Mamba/4.14.0-0

# If snakemake wants to regenerate outputs, do:
# rm -R .snakemake/ && snakemake --touch

# Create envs
# snakemake --use-conda --conda-frontend mamba --conda-create-envs-only

### Test
### Dry run
# snakemake -n --cores 1 --rerun-triggers mtime -f sim_output/MP/3.2.IQ_TREE/SEX~1.0/REC~1e-06/MUT~5e-07/{0..5}/sub/
# snakemake -n --cores 1 --rerun-triggers mtime --forcerun root_prune

### Real run
snakemake \
--executor slurm \
-j 80 \
--resources threads=${THREADS} mem_mb=$((5000 * 16)) \
--default-resources slurm_account=biol-bdelloids slurm_partition=devel \
--set-resources IQ_TREE:slurm_partition=short SINGER:slurm_partition=short \
--group-components subsample=40 iqtree=120 rootprune=80 singer=80 \
--use-conda --conda-frontend mamba \
--rerun-triggers mtime \
--rerun-incomplete \
--latency-wait 60000 \
--forcerun root_prune