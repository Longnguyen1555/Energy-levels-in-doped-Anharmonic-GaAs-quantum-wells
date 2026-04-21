function V = vconf_poly(z, V0_meV, beta1, beta2, k, e)
    if nargin < 6
        e = 1.602176634e-19;
    end
    V0_J = (V0_meV * 1e-3) * e;
    V = V0_J * (beta1 * (z / k).^2 + beta2 * (z / k).^8);
end
