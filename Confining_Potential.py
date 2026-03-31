import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

def vconf_poly(z, V0_meV, beta1, beta2, k, e=1.602176634e-19):

    V0_J = (V0_meV * 1e-3) * e
    return V0_J * (beta1 * (z / k) ** 2 + beta2 * (z / k) ** 8)