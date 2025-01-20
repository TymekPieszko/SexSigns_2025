import tskit, pyslim, msprime
import functions_smk as FS

ts_in = snakemake.input[0]
ts_out = snakemake.output[0]
vcf_out = snakemake.output[1]
MUT_RATE = float(snakemake.wildcards.MUT_RATE)
N_ANC = int(snakemake.params.N_ANC)
REC_RATE = 10 ** int(snakemake.params.REC_RATE)

ts = tskit.load(ts_in)
ts = pyslim.recapitate(ts, ancestral_Ne=N_ANC, recombination_rate=REC_RATE)
ts = FS.subsample_ts(ts, 2)
ts = msprime.sim_mutations(ts, rate=MUT_RATE)
ts.dump(ts_out)
with open(vcf_out, "w") as vcf_file:
    ts.write_vcf(vcf_file, allow_position_zero=True)
