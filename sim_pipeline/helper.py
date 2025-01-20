import re


def filter_targets(targets):
    new_targets = []
    for target in targets:
        SEX_FREQ = float(
            re.search(r"SEX_FREQ~([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", target).group(1)
        )
        GC_RATE = float(
            re.search(r"GC_RATE~([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", target).group(1)
        )
        if (
            SEX_FREQ != 0.0 and GC_RATE != 1.0e-6
        ):  # Want only one GC_RATE (1.0e-6) for non-zero SEX_FREQ.
            continue
        else:
            new_targets.append(target)
    return new_targets
