%% CODI PROJECTE ORBITAL
%---------------------------------------------
% Data: 23/05/2025
% Membres: Antonio Luque, Pau Cornudella i Alex Aleñà
% Assignatura: AERODINÀMICA, MECÀNICA DE VOL I ORBITAL
% Grup: 6
% ---------------------------------------------
clc; clear; close all;

%% 1. Leer datos.
[Year, DoY, Seconds, Constellation, SatID, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z, clock_offset] = leer_txt ('data.txt');

%Pasar de coordenadas cartesianas a geocentricas ECEF
[r, lambda_deg, phi_deg, lambda, phi_obs] = Cartesianes_geocentriques(x_TRF, y_TRF, z_TRF);

%% 2b 2c Posición satélite en funcion del tiempo
sat_list=r_vs_t(Seconds, SatID, r);

%2d) Velocidad satelite en funcion de la posición
v_vs_t(SatID, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z);

%2e) variació de rangs 
delta_r_km = delta_r(SatID, r);

%2f) comparacion 2015 vs 2025
[Year2025, DoY2025, Seconds2025, Constellation2025, SatID2025, x_TRF2025, y_TRF2025, z_TRF2025, v_x2025, v_y2025, v_z2025, clock_offset2025] = leer_txt ('data_2025.txt');
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

xlabel('Tiempo (h)');
ylabel('r (km)');
title('Comparación r(t) entre 2015 y 2025 (satélites comunes)');
legend show;

%% 3. Satellite ground tracks in the Terrestrial Reference Frame (TRF)
 satellite_groundtrack(lambda_deg, phi_deg, SatID);

%% 4. Constellation coverage in the Terrestrial Reference Frame (TRF)

unique_prns = unique(SatID,'stable');     % PRNs presentes
N_sats      = numel(unique_prns);         % nº de satélites
epochs      = unique(Seconds);            % vector de segundos únicos
N_epochs    = numel(epochs);              % nº de épocas (→ 288)

sat_data = cell(N_sats,1);
for k = 1:N_sats
    prn          = unique_prns(k);
    idx_prn      = (SatID == prn);
    [~,ord]      = sort(Seconds(idx_prn));  % ordenar por tiempo
    sat_data{k}.X = x_TRF(idx_prn); sat_data{k}.X = sat_data{k}.X(ord);
    sat_data{k}.Y = y_TRF(idx_prn); sat_data{k}.Y = sat_data{k}.Y(ord);
    sat_data{k}.Z = z_TRF(idx_prn); sat_data{k}.Z = sat_data{k}.Z(ord);
end

% Rejilla Fibonacci 
N_pts  = 1000;
grid   = Fibonacci(N_pts);          % puntos (x,y,z) unidad
lat_g  = asind(grid(:,3));                 % latitud  (°)
lon_g  = atan2d(grid(:,2), grid(:,1));     % longitud (°)
R_tierra = 6371000;                       % radio medio (m)

% Contadores por observador [ ≥1   ≥2   ≥3 ] satélites
vis_counts = zeros(N_pts,3);               % acumulado en épocas

for t = 1:N_epochs
    Xsat = zeros(N_sats,1); Ysat = Xsat; Zsat = Xsat;
    for s = 1:N_sats
        Xsat(s) = sat_data{s}.X(t);
        Ysat(s) = sat_data{s}.Y(t);
        Zsat(s) = sat_data{s}.Z(t);
    end

    for p = 1:N_pts
        obs = R_tierra * [ ...
              cosd(lat_g(p))*cosd(lon_g(p)); ...
              cosd(lat_g(p))*sind(lon_g(p)); ...
              sind(lat_g(p))               ];
        zen = obs / norm(obs);            
        n_vis = 0;

        for s = 1:N_sats
            los   = [Xsat(s);Ysat(s);Zsat(s)] - obs;
            elev  = asind( dot(los,zen)/norm(los) );   % elevación en °
            if elev > 5,  n_vis = n_vis + 1;  end
        end

        if n_vis>=1, vis_counts(p,1)=vis_counts(p,1)+1; end
        if n_vis>=2, vis_counts(p,2)=vis_counts(p,2)+1; end
        if n_vis>=3, vis_counts(p,3)=vis_counts(p,3)+1; end
    end
end

vis_hours = vis_counts * (5/60);           % N_pts × 3   (h)
[lat_grid, lon_grid] = meshgrid(linspace(-90,90,60), linspace(-180,180,120));
score = vis_hours(:,1) + 2*vis_hours(:,2) + 3*vis_hours(:,3);
vis_grid = griddata(lat_g, lon_g, score, lat_grid, lon_grid, 'linear');

figure;
contourf(lon_grid, lat_grid, vis_grid,20); shading interp
colormap(jet); colorbar;
xlabel('Longitud (°)'); ylabel('Latitud(°)');
title('Mapa de calor de la cobertura de Galileo ');
set(gca,'YDir','normal'); grid on;
xlim([-177 177]); ylim([-86  86]); %axis equal; grid on

[r_ECI, v_ECI, lambda_deg, phi_deg] = apartat_cinc_PAU(SatID, Seconds, DoY, Year, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z);