import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

from Helpers import *


def density_eq9(EF, E, Psi, p):

    kBT = p["kB"] * p["T"]
    pref = (p["mstar"] * kBT) / (np.pi * p["hbar"] ** 2)
    nj = pref * softplus((EF - E) / kBT)
    return (np.abs(Psi) ** 2) @ nj