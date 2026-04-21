function cache = precompute_thomas_eq20(n_interior)
% Precompute the constant coefficients of the tridiagonal matrix D
% from Cuesta 1995 Eq. (18), using the Thomas-style recurrence of Eq. (20).
% The matrix is: D = tridiag(1, -2, 1)

    if n_interior < 1
        error('n_interior must be >= 1');
    end

    a = ones(n_interior, 1);
    d = -2.0 * ones(n_interior, 1);
    c = ones(n_interior, 1);

    alpha = zeros(n_interior, 1);
    denom = zeros(n_interior, 1);

    denom(1) = d(1);
    if n_interior > 1
        alpha(1) = -c(1) / denom(1);
    else
        alpha(1) = 0.0;
    end

    for j = 2:n_interior
        denom(j) = d(j) + a(j) * alpha(j-1);
        if j < n_interior
            alpha(j) = -c(j) / denom(j);
        else
            alpha(j) = 0.0;
        end
    end

    cache.a = a;
    cache.d = d;
    cache.c = c;
    cache.alpha = alpha;
    cache.denom = denom;
end
