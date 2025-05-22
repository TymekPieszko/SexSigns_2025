import subprocess
import numpy as np


def subsample_ts(ts, ind_num):
    inds = np.random.choice(np.arange(ts.num_individuals), size=ind_num, replace=False)
    nodes = []
    for ind in inds:
        ind = ts.individual(ind)
        nodes.extend(ind.nodes)
    ts = ts.simplify(nodes, keep_input_roots=True)
    return ts, inds

def subsample_ts_to_given(ts, inds):
    nodes = []
    for ind in inds:
        ind = ts.individual(ind)
        nodes.extend(ind.nodes)
    ts = ts.simplify(nodes, keep_input_roots=True)
    return ts

# def get_windows(L, WINDOW_LEN):
#     windows = np.arange(0, L + WINDOW_LEN, WINDOW_LEN)
#     return windows


def get_sample_sets(IND_SAMPLE):
    sample_sets = [[i, i + 1] for i in range(0, 2 * IND_SAMPLE, 2)]
    return sample_sets


###############
### WITH LOG
###############
def run_singer(singer_path, Ne, m, vcf_prefix, arg_prefix, start, end, n, thin, log):
    cmd_args = [
        singer_path,
        "-Ne",
        Ne,
        "-m",
        m,
        "-vcf",
        vcf_prefix,
        "-output",
        arg_prefix,
        "-start",
        start,
        "-end",
        end,
        "-n",
        n,
        "-thin",
        thin,
    ]
    print(cmd_args)
    cmd = " ".join(cmd_args)
    subprocess.run(cmd, shell=True, stderr=log)


def run_convert_to_tskit(
    convert_to_tskit_path, arg_prefix, ts_prefix, start, end, step, log
):
    # convert_to_tskit -input prefix_of_arg_files -output prefix_of_tskit_files -start start_index -end end_index -step step_size
    cmd_args = [
        convert_to_tskit_path,
        "-input",
        arg_prefix,
        "-output",
        ts_prefix,
        "-start",
        start,
        "-end",
        end,
        "-step",
        step,
    ]
    cmd = " ".join(cmd_args)
    subprocess.run(cmd, shell=True, stderr=log)


# def run_singer(singer_path, Ne, m, vcf_prefix, arg_prefix, start, end, n, thin):
#     cmd_args = [
#         singer_path,
#         "-Ne",
#         Ne,
#         "-m",
#         m,
#         "-vcf",
#         vcf_prefix,
#         "-output",
#         arg_prefix,
#         "-start",
#         start,
#         "-end",
#         end,
#         "-n",
#         n,
#         "-thin",
#         thin,
#     ]
#     print(cmd_args)
#     cmd = " ".join(cmd_args)
#     subprocess.run(cmd, shell=True)


# def run_convert_to_tskit(
#     convert_to_tskit_path, arg_prefix, ts_prefix, start, end, step
# ):
#     # convert_to_tskit -input prefix_of_arg_files -output prefix_of_tskit_files -start start_index -end end_index -step step_size
#     cmd_args = [
#         convert_to_tskit_path,
#         "-input",
#         arg_prefix,
#         "-output",
#         ts_prefix,
#         "-start",
#         start,
#         "-end",
#         end,
#         "-step",
#         step,
#     ]
#     cmd = " ".join(cmd_args)
#     subprocess.run(cmd, shell=True)


def get_biased_intervals(ts, WINDOW_LEN, BIAS):
    windows = np.arange(0, ts.sequence_length + WINDOW_LEN, WINDOW_LEN)
    het_windows = ts.diversity(
        sample_sets=[[0, 1], [2, 3]], windows=windows, mode="site", span_normalise=True
    )
    het_windows = np.mean(
        het_windows, axis=1
    )  # Calc mean H_I per individual; take top BIAS proportion of windows
    index = int(len(het_windows) * BIAS)
    het_filter = np.argsort(het_windows)[-index:]
    intervals = [[i, j] for i, j in zip(windows[:-1], windows[1:])]
    intervals = np.array(intervals)[het_filter]
    return intervals


def infer_tree(fasta_f, tree_dir, threads):
    prefix = str(tree_dir) + "/" + fasta_f.stem.split(".")[0]
    cmd = (
        f"module load IQ-TREE/2.2.2.6-gompi-2022b && "
        f"iqtree2 -s {fasta_f} -m JC69 -pre {prefix} -T {threads} -redo -quiet"
    )
    subprocess.run(cmd, shell=True, check=True)
