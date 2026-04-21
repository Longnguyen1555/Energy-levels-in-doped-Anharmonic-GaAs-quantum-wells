function Psi_in = normalize_psi_interior(Psi_in, dz)
    norms = sqrt(sum(abs(Psi_in).^2, 1) * dz);
    Psi_in = Psi_in ./ norms;
end
