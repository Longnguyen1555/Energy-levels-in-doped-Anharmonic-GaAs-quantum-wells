% Main.m
% MATLAB version of the Python project

clear; clc;

Ldomain = 15e-9;
N = 1500;
z = linspace(-Ldomain/2.0, Ldomain/2.0, N).';

% Input data
Nd_sheet = 0;       % m^-2
doping_width = 2.0e-9; % 2 nm
B = 30;                 % Tesla
F = 0.0;               % V/m

V0_meV = 228.0;
beta1 = -2.0;
beta2 = 0.3;
k = 5e-9;
Vconf = vconf_poly(z, V0_meV, beta1, beta2, k);

% Physical parameters
p.e = 1.602176634e-19;
p.eps0 = 8.8541878128e-12;
m0 = 9.1093837015e-31;
p.epsr = 12.9;
p.mstar = 0.067 * m0;
p.hbar = 1.054571817e-34;
p.kB = 1.380649e-23;
p.T = 300.0;
p.B = B;
p.F = F;
p.L = Ldomain;
p.nStates = 4;
p.maxIter = 300;
p.mix = 0.5;          % Cuesta 1995 recommends lambda = 1/2 for stability.
p.tolEF = 1e-6 * p.e;
p.doping_width = doping_width;

out = sch_poisson_1d_cuesta1995(z, Nd_sheet, Vconf, p);

fprintf('Converged: %d\n', out.converged);
fprintf('Iterations: %d\n', out.iters);
fprintf('EF (meV): %.12g\n', out.EF / p.e * 1e3);
fprintf('Lowest energies (meV): ');
fprintf('%.12g ', out.E / p.e * 1e3);
fprintf('\n');

plot_levels_qw(z, Vconf, out.VH, out.E, out.Psi, p.e, true);
