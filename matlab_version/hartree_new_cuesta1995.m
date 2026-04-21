function VH = hartree_new_cuesta1995(z, dz, Nd_z, n_z, p, thomas_cache)
% Compute V_H^(new) using Cuesta 1995 section 3.2-3.3:
% - finite-difference Poisson equation -> D V_H = p
% - solve with the recurrence of Eqs. (19)-(20)

    rhs = (p.e^2 / (p.epsr * p.eps0)) * (Nd_z - n_z);
    p_vec = (dz^2) * rhs(2:end-1);

    n_interior = numel(p_vec);
    alpha = thomas_cache.alpha;
    denom = thomas_cache.denom;
    a = thomas_cache.a;

    gamma = zeros(n_interior, 1);
    gamma(1) = p_vec(1) / denom(1);

    for j = 2:n_interior
        gamma(j) = (p_vec(j) - a(j) * gamma(j-1)) / denom(j);
    end

    vh_in = zeros(n_interior, 1);
    vh_in(end) = gamma(end);

    for j = n_interior-1:-1:1
        vh_in(j) = alpha(j) * vh_in(j+1) + gamma(j);
    end

    VH = zeros(size(z));
    VH(2:end-1) = vh_in;
end
