from Equation_13 import *
from Confining_Potential import *
from Ploting import *

if __name__ == "__main__":
    # ====== CASE a) (Nd,B,F)=(0,0,0) ======
    Ldomain = 15e-9
    N = 1500
    z = np.linspace(-Ldomain/2, Ldomain/2, N)
    dz = z[1]-z[0]

    Nd_sheet = 1e17
    B = 0
    F = 0

    # Vconf example
    V0_meV = 228
    beta1, beta2 = -2.0, 0.3
    k = 5e-9
    Vconf = vconf_poly(z, V0_meV, beta1, beta2, k)

    # params
    e = 1.602176634e-19
    eps0 = 8.8541878128e-12
    m0 = 9.1093837015e-31
    p = dict(
        e=e, eps0=eps0, epsr=1,
        mstar=0.067*m0,
        hbar=1.054571817e-34,
        kB=1.380649e-23, T=4.2,
        B=B, F=F, L=Ldomain,
        nStates=4, maxIter=200, mix=0.3,
        tolEF=1e-6*e,
        useKineticPrefactor=True,
    )

    out = sch_poisson_1d(z, dz, Nd_sheet, Vconf, p)
    plot_levels_qw(z, Vconf, out["VH"], out["E"], out["Psi"], e_charge=e, plot_total_potential=False)