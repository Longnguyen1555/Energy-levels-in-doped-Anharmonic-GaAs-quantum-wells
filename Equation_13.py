import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq
from Equation_8 import *
from Equation_9 import *
from Equation_10 import *
from Helpers import *

def sch_poisson_1d_cuesta1995(z, Nd_sheet, Vconf, p):
    """
    Cuesta 1995 section 3.3.
    Step 7: compute V_H^(new) with the Eq. (19)-(20) recurrence.
    Step 9: mix V_H = lambda * V_H^(new) + (1-lambda) * V_H^(old).
    """
    dz = z[1] - z[0]
    N_full = z.size
    if N_full < 3:
        raise ValueError("Grid must contain at least 3 points including boundaries.")

    z_in = z[1:-1]
    Vconf_in = Vconf[1:-1]
    n_interior = z_in.size

    n_states = min(p["nStates"], max(1, n_interior - 2))
    mix = p["mix"]
    tolEF = p["tolEF"]

    Nd_z = N3d(z, Nd_sheet, p["doping_width"])
    thomas_cache = precompute_thomas_eq20(n_interior)

    VH_old = np.zeros(N_full)
    EF_old = 0.0
    converged = False

    for it in range(1, p["maxIter"] + 1):
        H = build_H_eq13(z_in, dz, Vconf_in, VH_old[1:-1], p)
        Psi_in, E = lowest_eigs(H, n_states)
        Psi_in = normalize_psi_interior(Psi_in, dz)
        Psi = pad_dirichlet_wavefunctions(Psi_in, N_full)

        if Nd_sheet <= 0.0:
            EF = np.min(E) - 100.0 * (p["kB"] * p["T"])
            n_z = np.zeros_like(z)
            VH_new = np.zeros_like(z)
        else:
            EF = solve_EF_eq10(Nd_sheet, E, p)
            n_z = density_eq9(EF, E, Psi, p)
            VH_new = hartree_new_cuesta1995(z, dz, Nd_z, n_z, p, thomas_cache)

        VH = mix * VH_new + (1.0 - mix) * VH_old

        if abs(EF - EF_old) < tolEF:
            converged = True
            VH_old = VH
            break

        EF_old = EF
        VH_old = VH

    return {
        "E": E,
        "Psi": Psi,
        "EF": EF,
        "n_z": n_z,
        "Nd_z": Nd_z,
        "VH": VH_old,
        "iters": it,
        "converged": converged,
    }
def build_H_eq13(z_in, dz, Vconf_in, VH_in, p):

    N = z_in.size
    alpha = p["hbar"] ** 2 / (2.0 * p["mstar"])

    off = -alpha / (dz ** 2) * np.ones(N - 1)
    Veff = (
        (p["e"] ** 2 * p["B"] ** 2) / (2.0 * p["mstar"]) * (z_in ** 2)
        - (p["e"] * p["F"]) * (z_in + p["L"] / 2.0)
        + Vconf_in
        + VH_in
    )
    main = 2.0 * alpha / (dz ** 2) + Veff
    return diags([off, main, off], offsets=[-1, 0, 1], shape=(N, N), format="csc")

def N3d(z, Nd_sheet, width):
    """
    Uniform N3d at z = 0 such that
    integral Nd(z) dz = Nd_sheet.
    """
    Nd_z = np.zeros_like(z)
    if Nd_sheet <= 0.0:
        return Nd_z

    if width <= 0.0:
        idx = np.argmin(np.abs(z))
        dz = z[1] - z[0]
        Nd_z[idx] = Nd_sheet / dz
        return Nd_z

    mask = np.abs(z) <= (width / 2.0)
    covered_width = np.trapezoid(mask.astype(float), z)
    if covered_width <= 0.0:
        idx = np.argmin(np.abs(z))
        dz = z[1] - z[0]
        Nd_z[idx] = Nd_sheet / dz
        return Nd_z

    Nd_z[mask] = Nd_sheet / covered_width
    return Nd_z