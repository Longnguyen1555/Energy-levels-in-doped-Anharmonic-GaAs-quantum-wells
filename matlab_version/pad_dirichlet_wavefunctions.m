function Psi = pad_dirichlet_wavefunctions(Psi_in, n_full)
    Psi = zeros(n_full, size(Psi_in, 2));
    Psi(2:end-1, :) = Psi_in;
end
