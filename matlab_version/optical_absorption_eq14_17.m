function [alpha12, alpha23, alpha13, hw_meV] = ...
    optical_absorption_eq14_17(E, Psi, EF, z, p)

%% CONSTANTS

mu0 = 4*pi*1e-7;
c   = 3e8;

nr = 3.2;
eps_opt = nr^2;

tau_in = 0.14e-12;

gamma = p.hbar / tau_in;

I_int = 0.5e10;

%% PHOTON ENERGY GRID

hw_meV = linspace(0,300,301);

hw_J = hw_meV * 1e-3 * p.e;

omega = hw_J / p.hbar;

%% DIPOLE MATRIX ELEMENTS

M = zeros(3);

for i = 1:3
    for j = 1:3

        M(i,j) = p.e * trapz( ...
            z, Psi(:,i).*z.*Psi(:,j));

    end
end

%% STORAGE

alpha12 = zeros(size(omega));
alpha23 = zeros(size(omega));
alpha13 = zeros(size(omega));

pairs = [1 2; 2 3; 1 3];

%% LOOP TRANSITIONS

for pp = 1:3

    i = pairs(pp,1);
    f = pairs(pp,2);

    %% ENERGY DIFFERENCE

    dE = E(f) - E(i);

    %% MATRIX ELEMENTS

    M_if = M(i,f);

    M_ii = M(i,i);

    M_ff = M(f,f);

    %% ==========================================================
    %% EQ.16
    %% sigma_if
    %% ==========================================================

    sigma_if = ...
    (p.mstar * p.kB * p.T) ...
    / ...
    (pi * p.hbar^2 * p.L) ...
    * ...
    log( ...
    (1 + exp((EF - E(i)) / (p.kB*p.T))) ...
    ./ ...
    (1 + exp((EF - E(f)) / (p.kB*p.T))) );

    %% ==========================================================
    %% EQ.14
    %% Linear absorption
    %% ==========================================================

    pre1 = omega .* ...
        sqrt(mu0 / (eps_opt * p.eps0));

    denom = (dE - hw_J).^2 + gamma^2;

    alpha1 = pre1 .* ...
        (abs(M_if)^2 .* sigma_if .* gamma) ...
        ./ denom;

    %% ==========================================================
    %% EQ.15
    %% Third-order nonlinear absorption
    %% ==========================================================

    pre3 = -2 .* omega .* ...
        sqrt(mu0 / (eps_opt * p.eps0)) .* ...
        (I_int / (p.eps0 * nr * c));

    bracket = 1 ...
        - (abs(M_ff - M_ii)^2 ...
        ./ (4 * abs(M_if)^2)) ...
        .* ...
        (((dE - hw_J).^2 ...
        - gamma^2 ...
        + 2*dE.*(dE - hw_J)) ...
        ./ ...
        (dE^2 + gamma^2));

    alpha3 = pre3 .* ...
        (abs(M_if)^4 .* sigma_if .* gamma) ...
        ./ ...
        (denom.^2) ...
        .* bracket;

    %% TOTAL ABSORPTION

    alpha_total = (alpha1 + alpha3) / 100;

    %% STORE

    if pp == 1

        alpha12 = alpha_total;

    elseif pp == 2

        alpha23 = alpha_total;

    else

        alpha13 = alpha_total;

    end
end

end