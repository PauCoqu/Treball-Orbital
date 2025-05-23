function delta_r_km = delta_r(SatID, r)


    sat_list = unique(SatID);
    delta_r_km = zeros(length(sat_list), 1);

    fprintf('Variación relativa de r (%%) para cada satélite:\n');
    for i = 1:length(sat_list)
        sat = sat_list(i);
        indices = SatID == sat;

        r_max = max(r(indices));
        r_min = min(r(indices));

        delta_r = (100 * (r_max - r_min) / r_min) / 1e3; 
        delta_r_km(i) = delta_r;

        fprintf('  PRN %d: %.3f %%\n', sat, delta_r);
    end
end


