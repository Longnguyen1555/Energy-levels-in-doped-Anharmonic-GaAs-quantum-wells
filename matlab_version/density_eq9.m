function n_z = density_eq9(EF, E, Psi, p)
    kBT = p.kB * p.T;
    pref = (p.mstar * kBT) / (pi * p.hbar^2);
    nj = pref * softplus((EF - E) / kBT);
    n_z = (abs(Psi).^2) * nj;
end
