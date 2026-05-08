function plot_figure3(F_array, E, e_charge, B_field, Nd_sheet)

    if nargin < 3 || isempty(e_charge)
        e_charge = 1.602176634e-19;
    end

    if nargin < 4
        B_field = 0;
    end

    if nargin < 5
        Nd_sheet = 0;
    end

    %% Convert units
    F_kVcm = F_array(:) / 1e5;

    % Joule -> meV
    E_meV = (E / e_charge) * 1e3;

    %% Check dimensions
    if size(E_meV,2) < 3
        error('E must contain at least 3 energy levels');
    end

    %% Energy separations
    dE21 = E_meV(:,2) - E_meV(:,1);
    dE31 = E_meV(:,3) - E_meV(:,1);
    dE32 = E_meV(:,3) - E_meV(:,2);

    %% Plot
    figure('Color','w');

    hold on;

    plot(F_kVcm, dE21, ...
        'b-', ...
        'LineWidth',2.5, ...
        'DisplayName','E_2 - E_1');

    plot(F_kVcm, dE31, ...
        'r:', ...
        'LineWidth',2.5, ...
        'DisplayName','E_3 - E_1');

    plot(F_kVcm, dE32, ...
        'k--', ...
        'LineWidth',2.5, ...
        'DisplayName','E_3 - E_2');

    xlabel('F (kV.cm^{-1})');
    ylabel('Energy separation (meV)');

    title('Energy separations vs electric field');

    grid on;
    box on;

    legend('Location','best');

    text(0.05,0.92, ...
        sprintf('B = %.1f T', B_field), ...
        'Units','normalized', ...
        'FontWeight','bold');

    text(0.05,0.86, ...
        sprintf('N_d = %.1e m^{-2}', Nd_sheet), ...
        'Units','normalized', ...
        'FontWeight','bold');

    hold off;

end