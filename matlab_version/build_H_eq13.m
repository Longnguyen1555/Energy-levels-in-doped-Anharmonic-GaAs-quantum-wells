function H = build_H_eq13(z_in, dz, Vconf_in, VH_in, p)
    N = numel(z_in);
    alpha = p.hbar^2 / (2.0 * p.mstar);

    z_in = z_in(:);
    Vconf_in = Vconf_in(:);
    VH_in = VH_in(:);

    off = (-alpha / dz^2) * ones(N-1, 1);
    Veff = (p.e^2 * p.B^2) / (2.0 * p.mstar) .* (z_in.^2) ...
         - (p.e * p.F) .* (z_in + p.L/2.0) ...
         + Vconf_in + VH_in;
    main = 2.0 * alpha / dz^2 + Veff;

    H = spdiags([[off;0], main, [0;off]], [-1,0,1], N, N);
    H = (H + H') / 2;
end
