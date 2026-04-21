function plot_levels_qw(z, Vconf, VH, E, Psi, e_charge, plot_total_potential)
    if nargin < 6 || isempty(e_charge)
        e_charge = 1.602176634e-19;
    end
    if nargin < 7
        plot_total_potential = true;
    end

    z_nm = z * 1e9;
    if plot_total_potential
        Vplot = Vconf + VH;
        lbl = 'Vconf + VH';
    else
        Vplot = Vconf;
        lbl = 'Vconf';
    end

    V_meV = (Vplot / e_charge) * 1e3;
    E_meV = (E / e_charge) * 1e3;

    n_show = min(4, numel(E_meV));
    Erange = max(V_meV) - min(V_meV);
    if Erange > 0.0
        amp = 0.15 * Erange;
    else
        amp = 10.0;
    end

    figure('Color', 'w');
    plot(z_nm, V_meV, 'LineWidth', 2, 'DisplayName', lbl);
    hold on;

    for i = 1:n_show
        prob = abs(Psi(:, i)).^2;
        if max(prob) > 0.0
            prob = prob / max(prob);
        end
        plot(z_nm, E_meV(i) * ones(size(z_nm)), '--', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('E%d', i));
        plot(z_nm, E_meV(i) + amp * prob, 'LineWidth', 2, ...
             'DisplayName', sprintf('|psi%d|^2', i));
    end

    xlabel('z (nm)');
    ylabel('Energy (meV)');
    grid on;
    title('Self-consistent levels with Cuesta 1995 Hartree update');
    legend('Location', 'best', 'FontSize', 9);
    hold off;
end
