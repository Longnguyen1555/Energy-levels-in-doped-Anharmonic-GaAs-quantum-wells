function EF = solve_EF_eq10(Nd_sheet, E, p)
%SOLVE_EF_EQ10  Solve the Fermi level from Eq. (10)
% using Halley's method (3rd-order Newton-Raphson).
%
%   EF = solve_EF_eq10(Nd_sheet, E, p)
%
%   Solves:
%
%       Nd_sheet = A * sum_i log(1 + exp(beta*(EF - E_i)))
%
%   where:
%
%       A    = mstar*kB*T/(pi*hbar^2)
%       beta = 1/(kB*T)
%
%   Energies are in Joule.

    % ============================================================
    % 1. Pre-processing
    % ============================================================

    E = real(E(:));

    if isempty(E)
        error('solve_EF_eq10:EmptyE', ...
              'Eigenvalue array E is empty.');
    end

    if any(~isfinite(E))
        error('solve_EF_eq10:NonFiniteE', ...
              'Eigenvalue array E contains NaN or Inf.');
    end

    % No doping case
    if Nd_sheet <= 0.0
        kBT = p.kB * p.T;
        EF = min(E) - 100.0 * kBT;
        return;
    end

    % ============================================================
    % 2. Constants
    % ============================================================

    kBT  = p.kB * p.T;

    A    = (p.mstar * kBT) / (pi * p.hbar^2);

    beta = 1.0 / kBT;

    % ============================================================
    % 3. Initial guess
    % ============================================================

    % Better initial guess than EF = 0
    EF = min(E);

    % ============================================================
    % 4. Tolerances
    % ============================================================

    if isfield(p, 'e') && isfinite(p.e) && p.e > 0.0
        tol_strict = 1e-12 * p.e;
    else
        tol_strict = 1e-12;
    end

    if isfield(p, 'tolEF') && isfinite(p.tolEF) && p.tolEF > 0.0
        tol = min(p.tolEF, tol_strict);
    else
        tol = tol_strict;
    end

    max_iter = 50;

    % ============================================================
    % 5. Halley iteration (3rd-order Newton-Raphson)
    % ============================================================

    for iter = 1:max_iter

        u = beta * (EF - E);

        % ---------- stable exp ----------
        u = min(u, 700.0);

        exp_u = exp(u);

        % ========================================================
        % f(EF)
        % ========================================================

        % f = A*sum(log(1+exp(u))) - Nd_sheet

        val_f = A * sum(stable_softplus(u)) - Nd_sheet;

        % ========================================================
        % f'(EF)
        % ========================================================

        % f' = A*beta*sum(sigmoid(u))

        sig = exp_u ./ (1.0 + exp_u);

        val_df = A * beta * sum(sig);

        % ========================================================
        % f''(EF)
        % ========================================================

        % f'' = A*beta^2*sum(exp(u)/(1+exp(u))^2)

        val_d2f = A * beta^2 * sum(exp_u ./ ((1.0 + exp_u).^2));

        % ========================================================
        % Safety checks
        % ========================================================

        if ~isfinite(val_f) || ...
           ~isfinite(val_df) || ...
           ~isfinite(val_d2f)

            error('solve_EF_eq10:NonFiniteHalley', ...
                  'Halley iteration produced NaN or Inf.');
        end

        if abs(val_df) < 1e-30
            error('solve_EF_eq10:SmallDerivative', ...
                  'First derivative too small.');
        end

        % ========================================================
        % Halley update
        % ========================================================

        denominator = 2.0 * (val_df^2) - val_f * val_d2f;

        if abs(denominator) < 1e-40
            error('solve_EF_eq10:SmallDenominator', ...
                  'Halley denominator too small.');
        end

        delta_EF = (2.0 * val_f * val_df) / denominator;

        % ========================================================
        % Step limiter
        % ========================================================

        max_step = 50.0 * kBT;

        if abs(delta_EF) > max_step
            delta_EF = sign(delta_EF) * max_step;
        end

        % Update
        EF = EF - delta_EF;

        % ========================================================
        % Convergence
        % ========================================================

        if abs(delta_EF) < tol
            return;
        end

        if abs(val_f) < 1e-12 * Nd_sheet
            return;
        end
    end

    % ============================================================
    % 6. Warning
    % ============================================================

    warning('solve_EF_eq10:NoConvergence', ...
            ['Halley method did not converge after %d iterations.'], ...
            max_iter);
end

% =================================================================
% Stable softplus
% =================================================================
function y = stable_softplus(x)

    y = max(x, 0.0) + log1p(exp(-abs(x)));

end