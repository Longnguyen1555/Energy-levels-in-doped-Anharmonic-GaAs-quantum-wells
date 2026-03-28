import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import brentq

def plot_levels_qw(z, Vconf, VH, E, Psi, e_charge=1.602176634e-19, plot_total_potential=False):
    z_nm = z * 1e9
    Vplot = Vconf + VH if plot_total_potential else Vconf
    V_meV = (Vplot / e_charge) * 1e3
    E_meV = (E     / e_charge) * 1e3

    nShow = min(4, len(E_meV))
    Erange = V_meV.max() - V_meV.min()
    amp = 0.15 * Erange

    plt.figure()
    plt.plot(z_nm, V_meV, linewidth=2, label="Vconf" if not plot_total_potential else "Vconf+VH")

    for i in range(nShow):
        prob = np.abs(Psi[:, i])**2
        prob = prob/prob.max()
        plt.plot(z_nm, E_meV[i]*np.ones_like(z_nm), "--", linewidth=1.2, label=f"E{i+1}")
        plt.plot(z_nm, E_meV[i] + amp*prob, linewidth=2, label=rf"|psi{i+1}|^2")

    plt.xlabel("z (nm)")
    plt.ylabel("Energy (meV)")
    plt.grid(True)
    plt.title("Energy levels, probabilities, and confining potential")
    plt.legend(loc="best")
    plt.show()