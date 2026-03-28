import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

def lowest_eigs(H, n_states):
    E, Psi = eigsh(H, k=n_states, which="SA")
    idx = np.argsort(E)
    return Psi[:, idx], E[idx]

def normalize_psi(Psi, dz):
    norms = np.sqrt(np.sum(np.abs(Psi)**2, axis=0) * dz)
    return Psi / norms

def softplus(x):
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)