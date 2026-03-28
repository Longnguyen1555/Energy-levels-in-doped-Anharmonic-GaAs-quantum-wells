import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

def hatree_potential(dz, Nd_z, n_z, p):
    rhs = (p["e"]**2/(p["epsr"]*p["eps0"])) * (Nd_z - n_z)
    N = rhs.size
    Ni = N - 2
    A = diags([np.ones(Ni), -2*np.ones(Ni), np.ones(Ni)], [-1, 0, 1], shape=(Ni, Ni), format="csc") / (dz**2)
    b = rhs[1:-1]
    VH = np.zeros(N)
    VH[1:-1] = spsolve(A, b)
    return VH
