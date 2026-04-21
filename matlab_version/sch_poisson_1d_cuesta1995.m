function out = sch_poisson_1d_cuesta1995(z, Nd_sheet, Vconf, p)
% Cuesta 1995 section 3.3.
% Step 7: compute V_H^(new) with the Eq. (19)-(20) recurrence.
% Step 9: mix V_H = lambda * V_H^(new) + (1-lambda) * V_H^(old).

    dz = z(2) - z(1);
    N_full = numel(z);
    if N_full < 3
        error('Grid must contain at least 3 points including boundaries.');
    end

    z_in = z(2:end-1);
    Vconf_in = Vconf(2:end-1);
    n_interior = numel(z_in);

    n_states = min(p.nStates, max(1, n_interior - 2));
    mix = p.mix;
    tolEF = p.tolEF;

    Nd_z = N3d(z, Nd_sheet, p.doping_width);
    thomas_cache = precompute_thomas_eq20(n_interior);

    VH_old = zeros(N_full, 1);
    EF_old = 0.0;
    converged = false;

    for it = 1:p.maxIter
        H = build_H_eq13(z_in, dz, Vconf_in, VH_old(2:end-1), p);
        [Psi_in, E] = lowest_eigs(H, n_states);
        Psi_in = normalize_psi_interior(Psi_in, dz);
        Psi = pad_dirichlet_wavefunctions(Psi_in, N_full);

        if Nd_sheet <= 0.0
            EF = min(E) - 100.0 * (p.kB * p.T);
            n_z = zeros(size(z));
            VH_new = zeros(size(z));
        else
            EF = solve_EF_eq10(Nd_sheet, E, p);
            n_z = density_eq9(EF, E, Psi, p);
            VH_new = hartree_new_cuesta1995(z, dz, Nd_z, n_z, p, thomas_cache);
        end

        VH = mix * VH_new + (1.0 - mix) * VH_old;

        if abs(EF - EF_old) < tolEF
            converged = true;
            VH_old = VH;
            break;
        end

        EF_old = EF;
        VH_old = VH;
    end

    out.E = E;
    out.Psi = Psi;
    out.EF = EF;
    out.n_z = n_z;
    out.Nd_z = Nd_z;
    out.VH = VH_old;
    out.iters = it;
    out.converged = converged;
end
