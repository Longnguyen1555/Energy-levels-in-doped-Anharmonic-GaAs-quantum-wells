function Nd_z = N3d(z, Nd_sheet, width)
% Uniform N3d at z = 0 such that integral Nd(z) dz = Nd_sheet.

    Nd_z = zeros(size(z));
    if Nd_sheet <= 0.0
        return;
    end

    if width <= 0.0
        [~, idx] = min(abs(z));
        dz = z(2) - z(1);
        Nd_z(idx) = Nd_sheet / dz;
        return;
    end

    mask = abs(z) <= (width / 2.0);
    covered_width = trapz(z, double(mask));
    if covered_width <= 0.0
        [~, idx] = min(abs(z));
        dz = z(2) - z(1);
        Nd_z(idx) = Nd_sheet / dz;
        return;
    end

    Nd_z(mask) = Nd_sheet / covered_width;
end
