function [Psi, E] = lowest_eigs(H, n_states)
    % Robust lowest eigenpairs for a real symmetric sparse Hamiltonian.
    H = (H + H') / 2;

    opts = struct();
    opts.issym = true;
    opts.isreal = true;
    opts.tol = 1e-10;
    opts.maxit = 2000;

    try
        [Psi, D] = eigs(H, n_states, 'smallestreal', opts);
    catch
        [Psi, D] = eigs(H, n_states, 0, opts);
    end

    E = real(diag(D));
    Psi = real(Psi);

    finite_mask = isfinite(E);
    if ~all(finite_mask)
        error('lowest_eigs:NonFiniteEigenvalues', ...
              'Non-finite eigenvalues returned by eigs.');
    end

    [E, idx] = sort(E, 'ascend');
    Psi = Psi(:, idx);
end
