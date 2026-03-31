from Equation_13 import *
from Confining_Potential import *
from Ploting import *

if __name__ == "__main__":

    Ldomain = 15e-9
    N = 1500
    z = np.linspace(-Ldomain / 2.0, Ldomain / 2.0, N)

    #Input data
    Nd_sheet = 1e17           # m^-2
    doping_width = 2.0e-9     # 2 nm
    B = 0                     # Tesla
    F = 0.0                   # V/m


    V0_meV = 228.0
    beta1, beta2 = -2.0, 0.3
    k = 5e-9
    Vconf = vconf_poly(z, V0_meV, beta1, beta2, k)

    # Physical parameters.
    e = 1.602176634e-19
    eps0 = 8.8541878128e-12
    m0 = 9.1093837015e-31
    p = dict(
        e=e,
        eps0=eps0,
        epsr=12.9,
        mstar=0.067 * m0,
        hbar=1.054571817e-34,
        kB=1.380649e-23,
        T=300.0,
        B=B,
        F=F,
        L=Ldomain,
        nStates=4,
        maxIter=300,
        mix=0.5,               # Cuesta 1995 recommends lambda = 1/2 for stability.
        tolEF=1e-6 * e,
        doping_width=doping_width,
    )

    out = sch_poisson_1d_cuesta1995(z, Nd_sheet, Vconf, p)

    print("Converged:", out["converged"])
    print("Iterations:", out["iters"])
    print("EF (meV):", out["EF"] / e * 1e3)
    print("Lowest energies (meV):", out["E"] / e * 1e3)

    plot_levels_qw(
        z,
        Vconf,
        out["VH"],
        out["E"],
        out["Psi"],
        e_charge=e,
        plot_total_potential=True,

    )