%% CODI PROJECTE ORBITAL
%---------------------------------------------
% Data: 23/05/2025
% Membres: Antonio Luque, Pau Cornudella i Alex Aleñà
% Assignatura: AERODINÀMICA, MECÀNICA DE VOL I ORBITAL
% Grup: 6
% ---------------------------------------------
clc; clear; close all;

[Year, DoY, Seconds, Constellation, SatID, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z, clock_offset] = llegir_txt ('data.txt');

[r, lambda_deg, phi_deg] = Cartesianes_geocentriques(x_TRF, y_TRF, z_TRF);

%2b 2c posicio satelit en funcio del temps
r_vs_t(Seconds, SatID, r);

%2d) velocitat satelit en funcio de la posició
v_vs_t(SatID, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z);

%2e) variació de rangs 
sat_list = unique(SatID);

for i = 1:length(sat_list)
    sat = sat_list(i);
    indices = SatID == sat;
    
    r_max = max(r(indices));
    r_min = min(r(indices));
    
    delta_r = (100*(r_max - r_min)/r_min) / 1e3; % en kilómetros
    delta_r_km(i) = delta_r;
    
    fprintf('  PRN %d: %.3f %\n', sat, delta_r);
end
