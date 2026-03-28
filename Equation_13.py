import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq
from Equation_8 import *
from Equation_9 import *
from Equation_10 import *
from Helpers import *

def sch_poisson_1d(z, dz, Nd_sheet, Vconf, p):
    N = z.size
    VH = np.zeros(N)
    EF = 0

    # donor profile for Poisson: uniform so integral Nd(z) dz = Nd_sheet
    Nd_z = np.full(N, Nd_sheet / (z.max() - z.min()))

    for it in range(1, p["maxIter"] + 1):
        EF_old = EF
        VH_old = VH.copy()

        H = build_matrix_Hij(z, dz, Vconf, VH_old, p)
        Psi, E = lowest_eigs(H, p["nStates"])
        Psi = normalize_psi(Psi, dz)

        if Nd_sheet <= 0:
            EF = E.min() - 100 * (p["kB"] * p["T"])
            n_z = np.zeros(N)
            VH = np.zeros(N)
            break

        EF = fermi_eq10(Nd_sheet, E, p)
        print(EF)
        n_z = density_eq9(EF, E, Psi, p)

        VH_new = hatree_potential(dz, Nd_z, n_z, p)
        VH = p["mix"] * VH_new + (1.0 - p["mix"]) * VH_old

        if abs(EF - EF_old) < p["tolEF"]:
            break

    return {"E": E, "Psi": Psi, "EF": EF, "n_z": n_z, "VH": VH, "iters": it}

def build_matrix_Hij(z, dz, Vconf, VH, p):
    N = z.size
    if p["useKineticPrefactor"]:
        alpha = p["hbar"]**2/(2*p["mstar"])
    else:
        alpha = 1.0

    off = -alpha/(dz**2) * np.ones(N)

    Veff = (p["e"]**2 * p["B"]**2)/(2*p["mstar"]) * (z**2) \
         - (p["e"]*p["F"]) * (z + p["L"]/2.0) \
         + Vconf + VH

    main = 2*alpha/(dz**2) + Veff
    return diags([off, main, off], offsets=[-1, 0, 1], shape=(N, N), format="csc")