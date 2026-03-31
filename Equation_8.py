import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

def hartree_new_cuesta1995(z, dz, Nd_z, n_z, p, thomas_cache):
    """
    Compute V_H^(new) using Cuesta 1995 section 3.2-3.3:

    - finite-difference Poisson equation -> D V_H = p
    - solve with the recurrence of Eqs. (19)-(20)
    """
    rhs = (p["e"] ** 2 / (p["epsr"] * p["eps0"])) * (Nd_z - n_z)
    p_vec = (dz ** 2) * rhs[1:-1]

    n_interior = p_vec.size
    alpha = thomas_cache["alpha"]
    denom = thomas_cache["denom"]
    a = thomas_cache["a"]

    gamma = np.zeros(n_interior)
    gamma[0] = p_vec[0] / denom[0]

    for j in range(1, n_interior):
        gamma[j] = (p_vec[j] - a[j] * gamma[j - 1]) / denom[j]

    vh_in = np.zeros(n_interior)
    vh_in[-1] = gamma[-1]

    for j in range(n_interior - 2, -1, -1):
        vh_in[j] = alpha[j] * vh_in[j + 1] + gamma[j]

    VH = np.zeros_like(z)
    VH[1:-1] = vh_in
    return VH

def precompute_thomas_eq20(n_interior):
    """
    Precompute the constant coefficients of the tridiagonal matrix D
    from Cuesta 1995 Eq. (18), using the Thomas-style recurrence of
    Eq. (20). The matrix is:

        D = tridiag(1, -2, 1)

    For the recurrence x_j = alpha_j * x_{j+1} + gamma_j,
    alpha depends only on D and can be precomputed once.
    """
    if n_interior < 1:
        raise ValueError("n_interior must be >= 1")

    a = np.ones(n_interior)
    d = -2.0 * np.ones(n_interior)
    c = np.ones(n_interior)

    alpha = np.zeros(n_interior)
    denom = np.zeros(n_interior)

    denom[0] = d[0]
    alpha[0] = (-c[0] / denom[0]) if n_interior > 1 else 0.0

    for j in range(1, n_interior):
        denom[j] = d[j] + a[j] * alpha[j - 1]
        alpha[j] = (-c[j] / denom[j]) if j < (n_interior - 1) else 0.0

    return {
        "a": a,
        "d": d,
        "c": c,
        "alpha": alpha,
        "denom": denom,
    }