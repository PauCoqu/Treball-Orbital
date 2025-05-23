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
delta_r_km = delta_r(SatID, r);




