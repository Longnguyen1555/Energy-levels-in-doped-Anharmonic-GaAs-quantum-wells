% Main.m
% MATLAB version of the Python project

clear; clc;

%% =======================================================================
% Choose figure
% ========================================================================

fprintf('=====================================\n');
fprintf(' Select figure to reproduce\n');
fprintf('=====================================\n');

fig_choice = input('Enter figure number: ');

%% Spatial grid
Ldomain = 15e-9;
N = 1500;

z = linspace(-Ldomain/2.0, Ldomain/2.0, N).';

%% Physical parameters
p.e = 1.602176634e-19;
p.eps0 = 8.8541878128e-12;

m0 = 9.1093837015e-31;

p.epsr = 12.9;
p.mstar = 0.067 * m0;

p.hbar = 1.054571817e-34;
p.kB = 1.380649e-23;
p.T = 300.0;

%% Structure parameters
V0_meV = 228.0;

beta1 = -2.0;
beta2 = 0.3;

k = 5e-9;

Vconf = vconf_poly(z, V0_meV, beta1, beta2, k);

%% Doping
Nd_sheet = 1e17;
doping_width = 2.0e-9;

%% External fields
B = 0;

%% Solver parameters
p.B = B;
p.F = 0;

p.L = Ldomain;

p.nStates = 4;

p.maxIter = 300;

p.mix = 0.5;
p.tolEF = 1e-6 * p.e;

p.doping_width = doping_width;

%% =======================================================================
% FIGURE 2
% =======================================================================

if fig_choice == 2

    fprintf('\nRunning Figure 2 simulation...\n');

    out = sch_poisson_1d_cuesta1995( ...
        z, Nd_sheet, Vconf, p);

    fprintf('Converged: %d\n', out.converged);
    fprintf('Iterations: %d\n', out.iters);

    fprintf('EF (meV): %.12g\n', ...
        out.EF / p.e * 1e3);

    fprintf('Lowest energies (meV): ');

    fprintf('%.12g ', ...
        out.E / p.e * 1e3);

    fprintf('\n');

    plot_levels_qw( ...
        z, ...
        Vconf, ...
        out.VH, ...
        out.E, ...
        out.Psi, ...
        p.e, ...
        true);

%% =======================================================================
% FIGURE 3
% =======================================================================

elseif fig_choice == 3

    fprintf('\nRunning Figure 3 simulation...\n');

    %% Electric-field sweep
    F_array = linspace(5e6, 25e6, 41);

    %% Storage
    E_results = zeros(length(F_array), p.nStates);

    for iF = 1:length(F_array)
    
        p.F = F_array(iF);
    
        fprintf('F = %.3g V/m (%d/%d)\n', ...
            p.F, iF, length(F_array));
    
        out = sch_poisson_1d_cuesta1995( ...
            z, Nd_sheet, Vconf, p);
    
        E_results(iF,:) = out.E(:).';
    
    end

    %% Plot Figure 3
    plot_figure3( ...
        F_array, ...
        E_results, ...
        p.e, ...
        p.B, ...
        Nd_sheet);

%% =======================================================================
% INVALID INPUT
% =======================================================================

elseif fig_choice == 5

    F_array = [50e5 150e5 250e5];

    alpha12_store = zeros(length(F_array),301);
    alpha23_store = zeros(length(F_array),301);
    alpha13_store = zeros(length(F_array),301);

    for s = 1:length(F_array)

        p.F = F_array(s);

        out = sch_poisson_1d_cuesta1995( ...
            z, Nd_sheet, Vconf, p);

        [a12, a23, a13, hw_meV] = ...
            optical_absorption_eq14_17( ...
            out.E, ...
            out.Psi, ...
            out.EF, ...
            z, ...
            p);

        alpha12_store(s,:) = a12;
        alpha23_store(s,:) = a23;
        alpha13_store(s,:) = a13;

    end

    plot_figure5( ...
        F_array, ...
        hw_meV, ...
        alpha12_store, ...
        alpha23_store, ...
        alpha13_store, ...
        p.B, ...
        Nd_sheet);

elseif fig_choice == 10

    fprintf('\n');
    fprintf('========================================\n');
    fprintf('   RECREATING FIGURE 10\n');
    fprintf('========================================\n');

    %% =====================================================
    %% MAGNETIC FIELD VALUES
    %% =====================================================

    B_array = [0 25 50];

    %% =====================================================
    %% STORAGE
    %% =====================================================

    alpha12_store = zeros(length(B_array),301);

    alpha23_store = zeros(length(B_array),301);

    alpha13_store = zeros(length(B_array),301);

    %% =====================================================
    %% LOOP OVER MAGNETIC FIELD
    %% =====================================================

    for s = 1:length(B_array)

        p.B = B_array(s);

        fprintf('Solving B = %.1f T ...\n', p.B);

        %% SELF-CONSISTENT SOLVER

        out = sch_poisson_1d_cuesta1995( ...
            z, ...
            Nd_sheet, ...
            Vconf, ...
            p);

        %% OPTICAL ABSORPTION

        [a12, a23, a13, hw_meV] = ...
            optical_absorption_eq14_17( ...
            out.E, ...
            out.Psi, ...
            out.EF, ...
            z, ...
            p);

        %% STORE

        alpha12_store(s,:) = a12;

        alpha23_store(s,:) = a23;

        alpha13_store(s,:) = a13;

    end

    %% =====================================================
    %% PLOT FIGURE 10
    %% =====================================================

    plot_figure10( ...
        B_array, ...
        hw_meV, ...
        alpha12_store, ...
        alpha23_store, ...
        alpha13_store, ...
        p.F, ...
        Nd_sheet);

    fprintf('\n');
    fprintf('Figure 10 completed successfully.\n'); 
    
else
  
    error('Invalid choice. Enter 2 or 3 .');

end