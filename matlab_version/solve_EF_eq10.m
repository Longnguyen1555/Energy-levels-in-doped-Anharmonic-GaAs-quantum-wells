function EF = solve_EF_eq10(Nd_sheet, E, p)
    kBT = p.kB * p.T;
    pref = (p.mstar * kBT) / (pi * p.hbar^2);

    E = real(E(:));
    if any(~isfinite(E))
        error('solve_EF_eq10:NonFiniteE', ...
              'Eigenvalue array E contains NaN or Inf.');
    end

    neutrality = @(EF) pref * sum(softplus((EF - E) / kBT)) - Nd_sheet;

    lo = min(E) - 50.0 * kBT;
    hi = max(E) + 50.0 * kBT;
    flo = neutrality(lo);
    fhi = neutrality(hi);

    cnt = 0;
    while cnt < 80 && (~isfinite(flo) || ~isfinite(fhi) || flo * fhi > 0.0)
        lo = lo - 50.0 * kBT;
        hi = hi + 50.0 * kBT;
        flo = neutrality(lo);
        fhi = neutrality(hi);
        cnt = cnt + 1;
    end

    if ~isfinite(lo) || ~isfinite(hi) || ~isfinite(flo) || ~isfinite(fhi)
        error('solve_EF_eq10:NonFiniteBracket', ...
              'Could not construct a finite bracket for fzero.');
    end

    if flo == 0.0
        EF = lo;
        return;
    elseif fhi == 0.0
        EF = hi;
        return;
    elseif flo * fhi > 0.0
        error('solve_EF_eq10:BracketFailure', ...
              'Could not bracket the Fermi level in solve_EF_eq10.');
    end

    EF = fzero(neutrality, [lo, hi]);
end
