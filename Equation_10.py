import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

from Helpers import *
def fermi_eq10(Nd_sheet, E, p):
    kBT = p["kB"] * p["T"]
    pref = (p["mstar"]*kBT)/(np.pi*p["hbar"]**2)

    def f(EF):
        return pref * np.sum(softplus((EF - E)/kBT)) - Nd_sheet

    lo = np.min(E) - 50*kBT
    hi = np.max(E) + 50*kBT
    flo, fhi = f(lo), f(hi)

    cnt = 0
    while flo*fhi > 0 and cnt < 80:
        lo -= 50*kBT
        hi += 50*kBT
        flo, fhi = f(lo), f(hi)
        cnt += 1

    return brentq(f, lo, hi)