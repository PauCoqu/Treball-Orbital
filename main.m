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

    % Crear figura
    figure;
    
    % Cargar imagen del mapa del mundo (ajusta el nombre si es distinto)
%     img = imread('World_map.png'); % Usa una proyección equirectangular
%     image([-180 180], [-90 90], flipud(img)); % Ajustar ejes a grados
%     set(gca,'YDir','normal'); % Invertir eje Y
%     hold on;
    
    % Ajustar ejes
    xlim([-180 180]);
    ylim([-90 90]);
    
    % Dibujar trayectorias
    plot(lambda_deg(SatID == 14), phi_deg(SatID == 14), 'r', 'LineWidth', 1.5, 'DisplayName', 'GSAT201');
    plot(lambda_deg(SatID == 18), phi_deg(SatID == 18), 'b', 'LineWidth', 1.5, 'DisplayName', 'GSAT202');
    
    xlabel('Longitud (°)');
    ylabel('Latitud (°)');
    title('Ground tracks de GSAT201 y GSAT202 (TRF)');
    legend show;

