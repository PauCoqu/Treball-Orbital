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

%2f) comparacio 2015 vs 2025
[Year2025, DoY2025, Seconds2025, Constellation2025, SatID2025, x_TRF2025, y_TRF2025, z_TRF2025, v_x2025, v_y2025, v_z2025, clock_offset2025] = llegir_txt ('data_2025.txt');
[r2025, lambda_deg2025, phi_deg2025] = Cartesianes_geocentriques(x_TRF2025, y_TRF2025, z_TRF2025);

r_vs_t(Seconds, SatID, r);
hold on

sat_comunes = intersect(unique(SatID), unique(SatID2025));

for i = 1:length(sat_comunes)
    sat = sat_comunes(i);
    indices = SatID2025 == sat;
    tiempo_horas_2025 = Seconds2025(indices) / 3600;
    r_km_2025 = r2025(indices) / 1e3;
    plot(tiempo_horas_2025, r_km_2025, '--', 'DisplayName', ['2025 - PRN ' num2str(sat)]);
end

xlabel('Temps (h)');
ylabel('r (km)');
title('Comparació r(t) entre 2015 y 2025 (satèl·lits comuns)');
legend show;

%3. Satellite ground tracks in the Terrestrial Reference Frame (TRF)

satellite_groundtrack(lambda_deg, phi_deg, SatID);

%4

% Paso 1: Crear vectores de latitudes y longitudes (cada 5 grados, por ejemplo)
lambda_4 = (-90:5:90)';      % vector columna
phi_4 = (-180:5:180)';   % vector columna

% Crear todos los pares posibles de (lon, lat)
[Lon4, Lat4] = meshgrid(lambda_4, phi_4);
Lon4 = Lon4(:);  % vector columna
Lat4 = Lat4(:);  % vector columna

% Ahora Lon y Lat son vectores (N x 1) donde N = número total de celdas

phi_rad4 = deg2rad(Lat4);       % latitud geocéntrica en radianes
lambda_rad4 = deg2rad(Lon4);    % longitud en radianes

x_obs = r(1) * cos(phi_rad4) .* cos(lambda_rad4);
y_obs = r(1) * cos(phi_rad4) .* sin(lambda_rad4);
z_obs = r(1) * sin(phi_rad4);

figure;
plot (phi_rad4, lambda_rad4,'.');

%Pasar de xyz a ENU

N = length(Lon4);  % número de observadores

[ENU]=xyz_2_ENU(x_obs, y_obs, z_obs, lambda_rad4, phi_rad4, N);

% [ENU_]
% 
% rho

plot3(ENU(1,:),ENU(2,:),ENU(3,:));



