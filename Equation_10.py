import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

from Helpers import *
def solve_EF_eq10(Nd_sheet, E, p):

    kBT = p["kB"] * p["T"]
    pref = (p["mstar"] * kBT) / (np.pi * p["hbar"] ** 2)

    def neutrality(EF):
        return pref * np.sum(softplus((EF - E) / kBT)) - Nd_sheet

    lo = np.min(E) - 50.0 * kBT
    hi = np.max(E) + 50.0 * kBT
    flo, fhi = neutrality(lo), neutrality(hi)

    cnt = 0
    while flo * fhi > 0.0 and cnt < 80:
        lo -= 50.0 * kBT
        hi += 50.0 * kBT
        flo, fhi = neutrality(lo), neutrality(hi)
        cnt += 1

    if flo * fhi > 0.0:
        raise RuntimeError("Could not bracket the Fermi level in solve_EF_eq10.")

    return brentq(neutrality, lo, hi)