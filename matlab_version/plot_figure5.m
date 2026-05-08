function plot_figure5( ...
    F_array, ...
    hw_meV, ...
    alpha12_store, ...
    alpha23_store, ...
    alpha13_store, ...
    B_field, ...
    Nd_sheet)

    %% Default arguments

    if nargin < 6
        B_field = 0;
    end

    if nargin < 7
        Nd_sheet = 0;
    end

    %% Convert field to kV/cm

    F_kVcm = F_array(:) / 1e5;

    %% Check dimensions

    if size(alpha12_store,1) ~= length(F_array)
        error('alpha12_store size mismatch');
    end

    if size(alpha23_store,1) ~= length(F_array)
        error('alpha23_store size mismatch');
    end

    if size(alpha13_store,1) ~= length(F_array)
        error('alpha13_store size mismatch');
    end

    %% Plot styles

    styles = {'-','--',':'};

    colors = {'b','k','r'};

    %% Labels

    labels = cell(length(F_array),1);

    for i = 1:length(F_array)

        labels{i} = sprintf( ...
            'F = %.0f kV cm^{-1}', ...
            F_kVcm(i));

    end

    %% =========================================================
    %% PLOT
    %% =========================================================

    figure('Color','w');

    %% ---------------------------------------------------------
    %% alpha12
    %% ---------------------------------------------------------

    subplot(3,1,1);

    hold on;
    box on;
    grid on;

    for s = 1:length(F_array)

        plot( ...
            hw_meV, ...
            alpha12_store(s,:), ...
            'Color', colors{1}, ...
            'LineStyle', styles{s}, ...
            'LineWidth', 2.5, ...
            'DisplayName', labels{s});

    end

    ylabel('\alpha_{12} (cm^{-1})');

    title('(a) \alpha_{12}');

    legend('Location','best');

    xlim([0 300]);

    %% ---------------------------------------------------------
    %% alpha23
    %% ---------------------------------------------------------

    subplot(3,1,2);

    hold on;
    box on;
    grid on;

    for s = 1:length(F_array)

        plot( ...
            hw_meV, ...
            alpha23_store(s,:), ...
            'Color', colors{2}, ...
            'LineStyle', styles{s}, ...
            'LineWidth', 2.5, ...
            'DisplayName', labels{s});

    end

    ylabel('\alpha_{23} (cm^{-1})');

    title('(b) \alpha_{23}');

    legend('Location','best');

    xlim([0 300]);

    %% ---------------------------------------------------------
    %% alpha13
    %% ---------------------------------------------------------

    subplot(3,1,3);

    hold on;
    box on;
    grid on;

    for s = 1:length(F_array)

        plot( ...
            hw_meV, ...
            alpha13_store(s,:), ...
            'Color', colors{3}, ...
            'LineStyle', styles{s}, ...
            'LineWidth', 2.5, ...
            'DisplayName', labels{s});

    end

    xlabel('Photon energy (meV)');

    ylabel('\alpha_{13} (cm^{-1})');

    title('(c) \alpha_{13}');

    legend('Location','best');

    xlim([0 300]);

    %% Global title

    sgtitle( ...
        'Variation of optical absorption coefficients');

    %% Extra information

    annotation( ...
        'textbox', ...
        [0.13 0.92 0.3 0.05], ...
        'String', sprintf('B = %.1f T', B_field), ...
        'EdgeColor', 'none', ...
        'FontWeight', 'bold');

    annotation( ...
        'textbox', ...
        [0.55 0.92 0.3 0.05], ...
        'String', sprintf('N_d = %.1e m^{-2}', Nd_sheet), ...
        'EdgeColor', 'none', ...
        'FontWeight', 'bold');

end