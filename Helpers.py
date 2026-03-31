import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

def lowest_eigs(H, n_states):
    E, Psi = eigsh(H, k=n_states, which="SA")
    idx = np.argsort(E)
    return Psi[:, idx], E[idx]

def normalize_psi_interior(Psi_in, dz):
    norms = np.sqrt(np.sum(np.abs(Psi_in) ** 2, axis=0) * dz)
    return Psi_in / norms

def softplus(x):
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)

def pad_dirichlet_wavefunctions(Psi_in, n_full):
    Psi = np.zeros((n_full, Psi_in.shape[1]), dtype=float)
    Psi[1:-1, :] = Psi_in
    return Psi

def softplus(x):
    """Stable log(1 + exp(x))."""
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)