# MATLAB version of Energy-levels-in-doped-Anharmonic-GaAs-quantum-wells

This MATLAB project is a line-by-line logical conversion of the original Python implementation.

## Files
- `Main.m`: entry point
- `sch_poisson_1d_cuesta1995.m`: self-consistent Schrodinger-Poisson loop
- `build_H_eq13.m`: Hamiltonian assembly
- `N3d.m`: donor density profile
- `hartree_new_cuesta1995.m`: Hartree update via Cuesta 1995 recurrence
- `precompute_thomas_eq20.m`: precompute Thomas-like coefficients
- `density_eq9.m`: electron density from occupied subbands
- `solve_EF_eq10.m`: Fermi level solver
- `vconf_poly.m`: confining potential
- `plot_levels_qw.m`: potential and level plotting
- `lowest_eigs.m`, `normalize_psi_interior.m`, `pad_dirichlet_wavefunctions.m`, `softplus.m`: helper utilities

## Run
Open MATLAB in this folder and run:

```matlab
Main
```

## Notes
- The logic and numerical flow are preserved from the Python version.
- Sparse matrices use `spdiags` and eigenpairs use `eigs`.
- Root finding uses `fzero` in place of SciPy `brentq`.
