import numpy as np
import tskit
import scipy.spatial.distance as ssd
import gzip
import subprocess
from Bio import SeqIO


def subsample_ts(ts, ind_num):
    inds = np.random.choice(ts.num_individuals, size=ind_num, replace=False)
    nodes = [n for i in inds for n in ts.individual(i).nodes]
    return ts.simplify(nodes)


def subsample_ts_to_given(ts, inds):
    nodes = []
    for ind in inds:
        ind = ts.individual(ind)
        nodes.extend(ind.nodes)
    ts = ts.simplify(nodes, keep_input_roots=True)
    return ts


def calc_hi(ts):
    hi_per_ind = ts.diversity(
        sample_sets=[i.nodes for i in ts.individuals()], mode="site"
    )
    hi = np.mean(hi_per_ind)
    return np.mean(hi)


def calc_hi_windows(ts, window_length):
    sample_sets = [i.nodes for i in ts.individuals()]
    num = int(ts.sequence_length / window_length) + 1
    windows = np.linspace(0, ts.sequence_length, num=num)
    hi_windows_per_ind = ts.diversity(
        sample_sets=sample_sets, windows=windows, mode="site"
    )
    hi_windows = np.mean(hi_windows_per_ind, axis=1)
    return hi_windows


def calc_fis(ts):
    fis_lst = []
    for var in ts.variants():
        gts = var.genotypes
        if set(gts) != {0, 1}:  # Only segregating, bi-allelic sites
            continue
        p = np.mean(gts)
        q = 1 - p
        hexp = 2 * p * q
        hobs = sum(a1 != a2 for a1, a2 in zip(gts[::2], gts[1::2])) / (len(gts) / 2)
        fis = 1 - hobs / hexp
        fis_lst.append(fis)
    return np.mean(fis_lst)


def calc_ld_ratio(ts, d1, d2):
    # Delete sites with > 1 mutation
    ids_to_remove = [site.id for site in ts.sites() if len(site.mutations) != 1]
    ts = ts.delete_sites(ids_to_remove)
    #####
    r2 = tskit.LdCalculator(ts).r2_matrix()
    positions = ts.sites_position
    distances = ssd.squareform(ssd.pdist(positions[:, None]))
    r2_d1 = np.nanmean(r2[distances < d1])
    r2_d2 = np.nanmean(r2[(distances > d1) & (distances < d2)])
    ld_ratio = r2_d1 / r2_d2
    return ld_ratio


def calc_r2(ts, d1, d2):
    # Delete sites with > 1 mutation
    ids_to_remove = [site.id for site in ts.sites() if len(site.mutations) != 1]
    ts = ts.delete_sites(ids_to_remove)
    #####
    r2 = tskit.LdCalculator(ts).r2_matrix()
    positions = ts.sites_position
    distances = ssd.squareform(ssd.pdist(positions[:, None]))
    r2 = np.nanmean(r2[(distances > d1) & (distances < d2)])
    return r2


def calc_r2_windows(ts, windows, target_sites):
    # Delete sites with > 1 mutation
    ids_to_remove = [site.id for site in ts.sites() if len(site.mutations) != 1]
    if target_sites > (ts.num_sites - len(ids_to_remove)):
        target_sites = ts.num_sites - len(ids_to_remove)
    ids_to_remove = np.concatenate(
        [
            ids_to_remove,
            np.random.choice(
                np.setdiff1d(np.arange(ts.num_sites), ids_to_remove),
                ts.num_sites - target_sites - len(ids_to_remove),
                replace=False,
            ),
        ]
    )
    ids_to_remove = [int(x) for x in ids_to_remove]
    # print("###############")
    print(ts.num_sites)
    # print(len(ids_to_remove))
    ts = ts.delete_sites(ids_to_remove)
    # print(ts.num_sites)
    # print("###############")
    r2 = tskit.LdCalculator(ts).r2_matrix()
    positions = ts.sites_position
    distances = ssd.squareform(ssd.pdist(positions[:, None]))
    r2_windows = [
        np.mean(r2[(distances > d1) & (distances < d2)])
        for d1, d2 in zip(windows[:-1], windows[1:])
    ]
    return r2_windows


def calc_delta_m(ts):
    dm_lst = []
    span_lst = []
    for tree in ts.trees():
        if tree.has_single_root:
            rank = tree.rank()
            if rank in [(4, 1), (4, 2)]:
                counts = [0, 0, 0, 0]
                for mut in tree.mutations():
                    if 0 <= mut.node <= 3:
                        counts[mut.node] += 1
                if rank == (4, 1):
                    w1 = abs(counts[0] - counts[2])
                    w2 = abs(counts[1] - counts[3])
                    b1 = abs(counts[0] - counts[1])
                    b2 = abs(counts[0] - counts[3])
                    b3 = abs(counts[2] - counts[1])
                    b4 = abs(counts[2] - counts[3])
                if rank == (4, 2):
                    w1 = abs(counts[0] - counts[3])
                    w2 = abs(counts[1] - counts[2])
                    b1 = abs(counts[0] - counts[1])
                    b2 = abs(counts[0] - counts[2])
                    b3 = abs(counts[3] - counts[1])
                    b4 = abs(counts[3] - counts[2])
                dm = ((b1 + b2 + b3 + b4) - 2 * (w1 + w2)) / np.sum(counts)
                dm_lst.append(dm)
                span_lst.append(tree.span)
    if len(dm_lst) == 0 or np.isnan(dm_lst).all():
        return np.nan
    mask = ~np.isnan(dm_lst)
    dm_lst = np.array(dm_lst)[mask]
    span_lst = np.array(span_lst)[mask]
    dm_ts = np.average(dm_lst, weights=span_lst)
    return dm_ts


def calc_delta_m_v2(ts):
    dm_trees = []
    for tree in ts.trees():
        if tree.has_single_root:
            rank = tree.rank()
            if rank in [(4, 1), (4, 2)]:
                counts = [0, 0, 0, 0]
                for mut in tree.mutations():
                    if 0 <= mut.node <= 3:
                        counts[mut.node] += 1
                if rank == (4, 1):
                    w1 = abs(counts[0] - counts[2])
                    w2 = abs(counts[1] - counts[3])
                    b1 = abs(counts[0] - counts[1])
                    b2 = abs(counts[0] - counts[3])
                    b3 = abs(counts[2] - counts[1])
                    b4 = abs(counts[2] - counts[3])
                if rank == (4, 2):
                    w1 = abs(counts[0] - counts[3])
                    w2 = abs(counts[1] - counts[2])
                    b1 = abs(counts[0] - counts[1])
                    b2 = abs(counts[0] - counts[2])
                    b3 = abs(counts[3] - counts[1])
                    b4 = abs(counts[3] - counts[2])
                dm = 0.25 * (b1 + b2 + b3 + b4) - 0.5 * (w1 + w2)
                dm_trees.append(dm)
    dm_ts = np.nanmean(dm_trees)
    return dm_ts


def calc_delta_m_iq(tree):
    desc = [
        [child.name for child in node.descendants if child.name is not None]
        for node in tree.walk()
    ]
    if len(desc[0]) == 0:  # Consider balanced trees
        lengths = {n.name: n.length for n in tree.walk() if n.name is not None}
        # print(lengths)
        tips = [set(x) for x in desc[1:] if len(x) == 2]
        # print(tips)
        if tips[0] == {"0", "2"} or tips[0] == {"1", "3"}:  # b0213
            w1 = abs(lengths["0"] - lengths["2"])
            w2 = abs(lengths["1"] - lengths["3"])
            b1 = abs(lengths["0"] - lengths["1"])
            b2 = abs(lengths["0"] - lengths["3"])
            b3 = abs(lengths["2"] - lengths["1"])
            b4 = abs(lengths["2"] - lengths["3"])
            dm = ((b1 + b2 + b3 + b4) - 2 * (w1 + w2)) / np.sum(list(lengths.values()))
            # if np.isnan(dm):
            #     print("Target tree gave NaN")
        elif tips[0] == {"0", "3"} or tips[0] == {"1", "2"}:  # b0312
            w1 = abs(lengths["0"] - lengths["3"])
            w2 = abs(lengths["1"] - lengths["2"])
            b1 = abs(lengths["0"] - lengths["1"])
            b2 = abs(lengths["0"] - lengths["2"])
            b3 = abs(lengths["3"] - lengths["1"])
            b4 = abs(lengths["3"] - lengths["2"])
            dm = ((b1 + b2 + b3 + b4) - 2 * (w1 + w2)) / np.sum(list(lengths.values()))
            # if np.isnan(dm):
            # print("Target tree gave NaN")
        else:  # Non-C balanced trees
            dm = np.nan
    else:
        dm = np.nan
    return dm


def get_tree_composition_21(ts, tree_dict):
    sample_sets = [[0, 1], [2, 3], [0, 2], [1, 3], [0, 3], [1, 2]]
    tracked = [tskit.Tree(ts, tracked_samples=samples) for samples in sample_sets]

    for i in range(ts.num_trees):
        for tree in tracked:
            tree.next()
            assert i == tree.index

        # Ignore trees with polytomies (retain trees with unary nodes)
        if any(
            tree.num_children(u) > 2
            for u in tracked[0].nodes()
            if tracked[0].is_internal(u)
        ):
            continue
        # Ignore trees with multiple roots
        if not tracked[0].has_single_root:
            continue

        nodes = list(tracked[0].nodes(order="timeasc"))
        assert nodes[-1] == tracked[0].root
        node4, node5 = nodes[4:6]

        counts = [
            [tree.num_tracked_samples(node4), tree.num_tracked_samples(node5)]
            for tree in tracked
        ]
        counts_set = [set(c) for c in counts]

        # BALANCED
        if counts_set[0] == counts_set[1] == {0, 2}:  # b0123
            if tracked[0].time(node4) == tracked[0].time(node5):
                tree_dict["S"]["b0123_eq"] += tracked[0].span
            else:
                if counts[0][0] == 2:
                    tree_dict["A"]["b0123_01"] += tracked[0].span
                if counts[0][0] == 0:
                    tree_dict["A"]["b0123_23"] += tracked[0].span
        if counts_set[2] == counts_set[3] == {0, 2}:  # b0213
            if tracked[0].time(node4) == tracked[0].time(node5):
                tree_dict["C"]["b0213_eq"] += tracked[0].span
            else:
                if counts[2][0] == 2:
                    tree_dict["S"]["b0213_02"] += tracked[0].span
                if counts[2][0] == 0:
                    tree_dict["S"]["b0213_13"] += tracked[0].span
        if counts_set[4] == counts_set[5] == {0, 2}:  # b0312
            if tracked[0].time(node4) == tracked[0].time(node5):
                tree_dict["C"]["b0312_eq"] += tracked[0].span
            else:
                if counts[4][0] == 2:
                    tree_dict["S"]["b0312_03"] += tracked[0].span
                if counts[4][0] == 0:
                    tree_dict["S"]["b0312_12"] += tracked[0].span
        # UNBALANCED
        if counts_set[0] == {2}:  # u01XX
            if counts_set[2] == counts_set[5] == {1, 2}:
                tree_dict["A"]["u0123"] += tracked[0].span
            if counts_set[3] == counts_set[4] == {1, 2}:
                tree_dict["A"]["u0132"] += tracked[0].span
        if counts_set[1] == {2}:  # u23XX
            if counts_set[2] == counts_set[4] == {1, 2}:
                tree_dict["A"]["u2301"] += tracked[0].span
            if counts_set[3] == counts_set[5] == {1, 2}:
                tree_dict["A"]["u2310"] += tracked[0].span
        if counts_set[2] == {2}:  # u02XX
            if counts_set[0] == counts_set[5] == {1, 2}:
                tree_dict["S"]["u0213"] += tracked[0].span
            if counts_set[1] == counts_set[4] == {1, 2}:
                tree_dict["S"]["u0231"] += tracked[0].span
        if counts_set[3] == {2}:  # u13XX
            if counts_set[0] == counts_set[4] == {1, 2}:
                tree_dict["S"]["u1302"] += tracked[0].span
            if counts_set[1] == counts_set[5] == {1, 2}:
                tree_dict["S"]["u1320"] += tracked[0].span
        if counts_set[4] == {2}:  # u03XX
            if counts_set[0] == counts_set[3] == {1, 2}:
                tree_dict["S"]["u0312"] += tracked[0].span
            if counts_set[1] == counts_set[2] == {1, 2}:
                tree_dict["S"]["u0321"] += tracked[0].span
        if counts_set[5] == {2}:  # u12XX
            if counts_set[0] == counts_set[2] == {1, 2}:
                tree_dict["S"]["u1203"] += tracked[0].span
            if counts_set[1] == counts_set[3] == {1, 2}:
                tree_dict["S"]["u1230"] += tracked[0].span
    return tree_dict


# def get_tree_composition_18(ts, tree_dict):
#     sample_sets = [[0, 1], [2, 3], [0, 2], [1, 3], [0, 3], [1, 2]]
#     tracked = [tskit.Tree(ts, tracked_samples=samples) for samples in sample_sets]

#     for i in range(ts.num_trees):
#         for tree in tracked:
#             tree.next()
#             assert i == tree.index

#         # Ignore trees with polytomies (retain trees with unary nodes)
#         if any(
#             tree.num_children(u) > 2
#             for u in tracked[0].nodes()
#             if tracked[0].is_internal(u)
#         ):
#             continue
#         # Ignore trees with multiple roots
#         if not tracked[0].has_single_root:
#             continue

#         nodes = list(tracked[0].nodes(order="timeasc"))
#         assert nodes[-1] == tracked[0].root
#         node4, node5 = nodes[4:6]

#         counts = [
#             [tree.num_tracked_samples(node4), tree.num_tracked_samples(node5)]
#             for tree in tracked
#         ]
#         counts_set = [set(c) for c in counts]

#         # BALANCED
#         if counts_set[0] == counts_set[1] == {0, 2}:  # b0123
#             if tracked[0].time(node4) == tracked[0].time(node5):
#                 pass
#             else:
#                 if counts[0][0] == 2:
#                     tree_dict["A"]["b0123_01"] += tracked[0].span
#                 if counts[0][0] == 0:
#                     tree_dict["A"]["b0123_23"] += tracked[0].span
#         if counts_set[2] == counts_set[3] == {0, 2}:  # b0213
#             if tracked[0].time(node4) == tracked[0].time(node5):
#                 pass
#             else:
#                 if counts[2][0] == 2:
#                     tree_dict["C"]["b0213_02"] += tracked[0].span
#                 if counts[2][0] == 0:
#                     tree_dict["C"]["b0213_13"] += tracked[0].span
#         if counts_set[4] == counts_set[5] == {0, 2}:  # b0312
#             if tracked[0].time(node4) == tracked[0].time(node5):
#                 pass
#             else:
#                 if counts[4][0] == 2:
#                     tree_dict["C"]["b0312_03"] += tracked[0].span
#                 if counts[4][0] == 0:
#                     tree_dict["C"]["b0312_12"] += tracked[0].span
#         # UNBALANCED
#         if counts_set[0] == {2}:  # u01XX
#             if counts_set[2] == counts_set[5] == {1, 2}:
#                 tree_dict["A"]["u0123"] += tracked[0].span
#             if counts_set[3] == counts_set[4] == {1, 2}:
#                 tree_dict["A"]["u0132"] += tracked[0].span
#         if counts_set[1] == {2}:  # u23XX
#             if counts_set[2] == counts_set[4] == {1, 2}:
#                 tree_dict["A"]["u2301"] += tracked[0].span
#             if counts_set[3] == counts_set[5] == {1, 2}:
#                 tree_dict["A"]["u2310"] += tracked[0].span
#         if counts_set[2] == {2}:  # u02XX
#             if counts_set[0] == counts_set[5] == {1, 2}:
#                 tree_dict["S"]["u0213"] += tracked[0].span
#             if counts_set[1] == counts_set[4] == {1, 2}:
#                 tree_dict["S"]["u0231"] += tracked[0].span
#         if counts_set[3] == {2}:  # u13XX
#             if counts_set[0] == counts_set[4] == {1, 2}:
#                 tree_dict["S"]["u1302"] += tracked[0].span
#             if counts_set[1] == counts_set[5] == {1, 2}:
#                 tree_dict["S"]["u1320"] += tracked[0].span
#         if counts_set[4] == {2}:  # u03XX
#             if counts_set[0] == counts_set[3] == {1, 2}:
#                 tree_dict["S"]["u0312"] += tracked[0].span
#             if counts_set[1] == counts_set[2] == {1, 2}:
#                 tree_dict["S"]["u0321"] += tracked[0].span
#         if counts_set[5] == {2}:  # u12XX
#             if counts_set[0] == counts_set[2] == {1, 2}:
#                 tree_dict["S"]["u1203"] += tracked[0].span
#             if counts_set[1] == counts_set[3] == {1, 2}:
#                 tree_dict["S"]["u1230"] += tracked[0].span
#     return tree_dict


def get_tree_dict_21():
    tree_dict = {
        "C": {
            "b0213_eq": 0,
            "b0312_eq": 0,
        },
        "A": {
            "u0123": 0,
            "u0132": 0,
            "u2301": 0,
            "u2310": 0,
            "b0123_01": 0,
            "b0123_23": 0,
        },
        "S": {
            "b0123_eq": 0,
            "b0213_02": 0,
            "b0213_13": 0,
            "b0312_03": 0,
            "b0312_12": 0,
            "u0213": 0,
            "u0312": 0,
            "u1203": 0,
            "u1302": 0,
            "u0231": 0,
            "u0321": 0,
            "u1230": 0,
            "u1320": 0,
        },
    }
    return tree_dict


# def get_tree_dict_18():
#     tree_dict = {
#         "C": {
#             "b0213_02": 0,
#             "b0213_13": 0,
#             "b0312_03": 0,
#             "b0312_12": 0,
#         },
#         "A": {
#             "u0123": 0,
#             "u0132": 0,
#             "u2301": 0,
#             "u2310": 0,
#             "b0123_01": 0,
#             "b0123_23": 0,
#         },
#         "S": {
#             "u0213": 0,
#             "u0312": 0,
#             "u1203": 0,
#             "u1302": 0,
#             "u0231": 0,
#             "u0321": 0,
#             "u1230": 0,
#             "u1320": 0,
#         },
#     }
#     return tree_dict


def get_biased_intervals(ts, window, bias):
    windows = np.arange(0, ts.sequence_length + window, window)
    het = ts.diversity(
        [[0, 1], [2, 3]], windows=windows, mode="site", span_normalise=True
    )
    het = np.mean(het, axis=1)
    intervals = np.column_stack([windows[:-1], windows[1:]])
    num_windows = int(len(het) * bias)
    indices = np.sort(np.argsort(het)[-num_windows:])
    return intervals[indices]


def rename_nodes(tree):
    # This can be simplied using a newick method rename
    node_names = [node.name for node in tree.walk() if node.name is not None]

    def get_number(name):
        return int(name.replace("n", ""))

    node_names = sorted(node_names, key=get_number)
    name_dict = {name: str(i) for i, name in enumerate(node_names)}
    # print(name_dict)
    for node in tree.walk():
        if node.name is not None:
            node.name = name_dict[node.name]
    return tree


def get_tree_category(tree, window):
    tree = rename_nodes(tree)
    # print(tree.ascii_art())
    desc = [
        [child.name for child in node.descendants if child.name is not None]
        for node in tree.walk()
    ]
    # print(desc)
    if len(desc[0]) == 1:
        category = "unbalanced"
    elif len(desc[0]) == 0:
        category = "balanced"
    # print(category)
    if category == "unbalanced":
        if desc[0][0] == "0":
            test_node = [x for x in desc[1:] if len(x) == 1][0][0]
            if test_node == "1":
                cat_dict = {"A": {"u2310": window}}
            elif test_node == "2":
                cat_dict = {"S": {"u1320": window}}
            elif test_node == "3":
                cat_dict = {"S": {"u1230": window}}
        elif desc[0][0] == "1":
            test_node = [x for x in desc[1:] if len(x) == 1][0][0]
            if test_node == "0":
                cat_dict = {"A": {"u2301": window}}
            elif test_node == "2":
                cat_dict = {"S": {"u0321": window}}
            elif test_node == "3":
                cat_dict = {"S": {"u0231": window}}
        elif desc[0][0] == "2":
            test_node = [x for x in desc[1:] if len(x) == 1][0][0]
            if test_node == "0":
                cat_dict = {"S": {"u1302": window}}
            elif test_node == "1":
                cat_dict = {"S": {"u0312": window}}
            elif test_node == "3":
                cat_dict = {"A": {"u0132": window}}
        if desc[0][0] == "3":
            test_node = [x for x in desc[1:] if len(x) == 1][0][0]
            if test_node == "0":
                cat_dict = {"S": {"u1203": window}}
            elif test_node == "1":
                cat_dict = {"S": {"u0213": window}}
            elif test_node == "2":
                cat_dict = {"A": {"u0123": window}}
    elif category == "balanced":
        test_nodes = [set(x) for x in desc[1:] if len(x) == 2]
        test_times = {n.name: n.length for n in tree.walk() if n.name is not None}
        # print(test_nodes)
        if test_nodes[0] == {"0", "1"} or test_nodes[0] == {"2", "3"}:
            t01 = test_times["0"] + test_times["1"]
            t23 = test_times["2"] + test_times["3"]
            if t01 < t23:
                cat_dict = {"A": {"b0123_01": window}}
            elif t01 > t23:
                cat_dict = {"A": {"b0123_23": window}}
            elif t01 == t23:
                cat_dict = np.random.choice(
                    [{"A": {"b0123_01": window}}, {"A": {"b0123_23": window}}]
                )
        elif test_nodes[0] == {"0", "2"} or test_nodes[0] == {"1", "3"}:
            t02 = test_times["0"] + test_times["2"]
            t13 = test_times["1"] + test_times["3"]
            if t02 < t13:
                cat_dict = {"S": {"b0213_02": window}}
            elif t02 > t13:
                cat_dict = {"S": {"b0213_13": window}}
            elif t02 == t13:
                cat_dict = np.random.choice(
                    [{"S": {"b0213_02": window}}, {"S": {"b0213_13": window}}]
                )
        elif test_nodes[0] == {"0", "3"} or test_nodes[0] == {"1", "2"}:
            t03 = test_times["0"] + test_times["3"]
            t12 = test_times["1"] + test_times["2"]
            if t03 < t12:
                cat_dict = {"S": {"b0312_03": window}}
            elif t03 > t12:
                cat_dict = {"S": {"b0312_12": window}}
            elif t03 == t12:
                cat_dict = np.random.choice(
                    [{"S": {"b0312_03": window}}, {"S": {"b0312_12": window}}]
                )
    return cat_dict


def calc_hi_from_fasta(fasta_f, samples):
    def get_number(sample):
        return int(sample[1:])

    samples.sort(key=get_number)
    with gzip.open(fasta_f, "rt") as f:
        fasta = list(SeqIO.parse(f, "fasta"))
    seqs = [rec.seq for rec in fasta if rec.name in samples]
    diff1 = sum(a != b for a, b in zip(seqs[0], seqs[1]))
    diff2 = sum(a != b for a, b in zip(seqs[2], seqs[3]))
    diff1, diff2 = diff1 / len(seqs[0]), diff2 / len(seqs[0])
    hi = (diff1 + diff2) / 2
    return hi


def filter_by_hi(cat_dicts, hi_windows, bias):
    index = int(len(hi_windows) * bias)
    hi_filter = np.argsort(hi_windows)[-index:]
    cat_dicts = [cat_dicts[i] for i in hi_filter]
    return cat_dicts


def sum_cat_dicts(cat_dicts):
    tree_dict = get_tree_dict_21()
    for cat_dict in cat_dicts:
        for k1, v1 in cat_dict.items():
            for k2, v2 in v1.items():
                tree_dict[k1][k2] += v2
    return tree_dict


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


def calc_class_props(tree_dict):
    C = (
        sum(tree_dict["C"].values())
        + tree_dict["S"]["b0213_02"]
        + tree_dict["S"]["b0213_13"]
        + tree_dict["S"]["b0312_03"]
        + tree_dict["S"]["b0312_12"]
    )
    S = (
        sum(tree_dict["S"].values())
        - tree_dict["S"]["b0213_02"]
        - tree_dict["S"]["b0213_13"]
        - tree_dict["S"]["b0312_03"]
        - tree_dict["S"]["b0312_12"]
        - tree_dict["S"]["b0123_eq"]
    )
    A = sum(tree_dict["A"].values()) + tree_dict["S"]["b0123_eq"]
    total = C + A + S
    if total == 0:
        return {"C": np.nan, "A": np.nan, "S": np.nan}
    return {"C": C / total, "A": A / total, "S": S / total}
