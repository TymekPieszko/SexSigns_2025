import numpy as np


def calc_hi(MUT, gamma):
    hi = (2 * MUT) / ((2 * MUT) + gamma)
    return hi


def calc_gamma_MP(REC, TRACT):
    gamma = REC * TRACT
    return gamma


def calc_gamma_CF_CI(REC, pos):
    gamma = (1 / 4) * (1 - np.exp(-2 * pos * REC))
    return gamma


def calc_gamma_CF_NCI(REC, pos):
    gamma = (1 / 3) * (1 - np.exp(-1.5 * pos * REC))
    return gamma