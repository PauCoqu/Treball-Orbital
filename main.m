%% CODI PROJECTE ORBITAL
%---------------------------------------------
% Data: 23/05/2025
% Membres: Antonio Luque, Pau Cornudella i Alex Aleñà
% Assignatura: AERODINÀMICA, MECÀNICA DE VOL I ORBITAL
% Grup: 6
% ---------------------------------------------
clc; clear; close all;

%1. Leer datos.
[Year, DoY, Seconds, Constellation, SatID, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z, clock_offset] = leer_txt ('data.txt');

[r, lambda_deg, phi_deg, lambda, phi_obs] = Cartesianes_geocentriques(x_TRF, y_TRF, z_TRF);

%2b 2c posicio satelit en funcio del temps
sat_list=r_vs_t(Seconds, SatID, r);

%2d) velocitat satelit en funcio de la posició
v_vs_t(SatID, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z);

%2e) variació de rangs 
delta_r_km = delta_r(SatID, r);

%2f) comparacio 2015 vs 2025
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

% %3. Satellite ground tracks in the Terrestrial Reference Frame (TRF)
% 
 satellite_groundtrack(lambda_deg, phi_deg, SatID);
% 
% %4
% 
% %% ------------------------------------------------------------------------
% %% 4. Constellation coverage in the Terrestrial Reference Frame (TRF)
% disp('Question 4: Assessing constellation coverage with a Fibonacci lattice');
% 
% % 4-A  Agrupar muestras por PRN y por época --------------------------------
% unique_prns = unique(SatID,'stable');     % PRNs presentes
% N_sats      = numel(unique_prns);         % nº de satélites
% epochs      = unique(Seconds);            % vector de segundos únicos
% N_epochs    = numel(epochs);              % nº de épocas (→ 288)
% 
% % Celda para guardar XYZ (filas = época)
% sat_data = cell(N_sats,1);
% for k = 1:N_sats
%     prn          = unique_prns(k);
%     idx_prn      = (SatID == prn);
%     [~,ord]      = sort(Seconds(idx_prn));  % ordenar por tiempo
%     sat_data{k}.X = x_TRF(idx_prn); sat_data{k}.X = sat_data{k}.X(ord);
%     sat_data{k}.Y = y_TRF(idx_prn); sat_data{k}.Y = sat_data{k}.Y(ord);
%     sat_data{k}.Z = z_TRF(idx_prn); sat_data{k}.Z = sat_data{k}.Z(ord);
% end
% 
% % 4-B  Rejilla Fibonacci (1 000 observadores distribuidos casi uniformes)
% N_pts  = 1000;
% grid   = fibonacci_sphere(N_pts);          % puntos (x,y,z) unidad
% lat_g  = asind(grid(:,3));                 % latitud  (°)
% lon_g  = atan2d(grid(:,2), grid(:,1));     % longitud (°)
% R_earth = 6371000;                       % radio medio (m)
% 
% % 4-C  Contadores por observador [ ≥1   ≥2   ≥3 ] satélites
% vis_counts = zeros(N_pts,3);               % acumulado en épocas
% 
% % 4-D  Bucle por época y por punto de la rejilla ---------------------------
% for t = 1:N_epochs
%     % XYZ de los N_sats en la época t
%     Xsat = zeros(N_sats,1); Ysat = Xsat; Zsat = Xsat;
%     for s = 1:N_sats
%         Xsat(s) = sat_data{s}.X(t);
%         Ysat(s) = sat_data{s}.Y(t);
%         Zsat(s) = sat_data{s}.Z(t);
%     end
% 
%     % Bucle por cada observador
%     for p = 1:N_pts
%         obs = R_earth * [ ...
%               cosd(lat_g(p))*cosd(lon_g(p)); ...
%               cosd(lat_g(p))*sind(lon_g(p)); ...
%               sind(lat_g(p))               ];
%         zen = obs / norm(obs);             % vector Up local
%         n_vis = 0;
% 
%         for s = 1:N_sats
%             los   = [Xsat(s);Ysat(s);Zsat(s)] - obs;
%             elev  = asind( dot(los,zen)/norm(los) );   % elevación en °
%             if elev > 5,  n_vis = n_vis + 1;  end
%         end
% 
%         % Acumular: ≥1, ≥2, ≥3
%         if n_vis>=1, vis_counts(p,1)=vis_counts(p,1)+1; end
%         if n_vis>=2, vis_counts(p,2)=vis_counts(p,2)+1; end
%         if n_vis>=3, vis_counts(p,3)=vis_counts(p,3)+1; end
%     end
% end
% 
% % 4-E  Pasar épocas (5 min) a horas
% vis_hours = vis_counts * (5/60);           % N_pts × 3   (h)
% 
% % 4-F  Interpolar sobre una retícula regular para el mapa
% [lat_grid, lon_grid] = meshgrid(linspace(-90,90,60), linspace(-180,180,120));
% score = vis_hours(:,1) + 2*vis_hours(:,2) + 3*vis_hours(:,3);
% vis_grid = griddata(lat_g, lon_g, score, lat_grid, lon_grid, 'linear');
% 
% figure;
% pcolor(lon_grid, lat_grid, vis_grid); shading interp
% colormap(jet); colorbar;
% xlabel('Longitude (°)'); ylabel('Latitude (°)');
% title('Coverage heat-map (ponderada por nº de satélites · h)');
% set(gca,'YDir','normal'); grid on;
% 
% %% ---------- ya debes tener la figura con la esfera y los puntos ----------
% figure, hold on, axis equal off
% view(40,25), colormap(jet)
% lighting gouraud, camlight headlight
% 
% % -------- esfera base ------------
% [xe,ye,ze] = sphere(120);
% surf(R_earth*xe, R_earth*ye, R_earth*ze, ...
%      'FaceColor',[0.4 0.7 0.9], 'EdgeColor','none', 'FaceAlpha',0.25);
% 
% % ---------- convertir lat/lon de la malla Fibonacci a XYZ --------------
% xg = R_earth * cosd(lat_g) .* cosd(lon_g);
% yg = R_earth * cosd(lat_g) .* sind(lon_g);
% zg = R_earth * sind(lat_g);
% 
% % -------- puntos coverage --------
% scatter3(xg, yg, zg, 28, score, 'filled')
% caxis([min(score) max(score)]), cb = colorbar;
% cb.Label.String = 'Horas ponderadas';
% title('Cobertura Galileo (>5°) con contorno mundial y ejes')
% 
% %% ==================  Ejes cartesianos ==================
% L = 1.1*R_earth;                    % longitud de las flechas
% plot3([0 L],[0 0],[0 0],'k','LineWidth',1.6)          % eje X
% plot3([0 0],[0 L],[0 0],'k','LineWidth',1.6)          % eje Y
% plot3([0 0],[0 0],[0 L],'k','LineWidth',1.6)          % eje Z
% text(L,0,0,'  X','FontWeight','bold')
% text(0,L,0,'  Y','FontWeight','bold')
% text(0,0,L,'  Z','FontWeight','bold')
% 
% %% ==================  Contorno del mundo =================
% load coastlines              % coastlat, coastlon (en grados)
% % Convertir a radianes y a XYZ
% clat = deg2rad(coastlat);
% clon = deg2rad(coastlon);
% x_coast = 1.002*R_earth*cos(clat).*cos(clon);  % 1.002 para que sobresalga
% y_coast = 1.002*R_earth*cos(clat).*sin(clon);
% z_coast = 1.002*R_earth*sin(clat);
% 
% plot3(x_coast, y_coast, z_coast, 'k', 'LineWidth',0.7)
% 
% %% Opcional: cuadrícula de meridianos / paralelos finos
% for latLine = -60:30:60
%     latR = deg2rad(latLine);
%     lon = linspace(-180,180,360);
%     xq = 1.001*R_earth*cosd(latLine).*cosd(lon);
%     yq = 1.001*R_earth*cosd(latLine).*sind(lon);
%     zq = 1.001*R_earth*sind(latLine)*ones(size(lon));
%     plot3(xq,yq,zq,'Color',[.5 .5 .5],'LineStyle',':')
% end
% for lonLine = -150:30:150
%     lonR = deg2rad(lonLine);
%     lat = linspace(-90,90,180);
%     xq = 1.001*R_earth*cosd(lat).*cosd(lonLine);
%     yq = 1.001*R_earth*cosd(lat).*sind(lonLine);
%     zq = 1.001*R_earth*sind(lat);
%     plot3(xq,yq,zq,'Color',[.5 .5 .5],'LineStyle',':')
% end
% 
% rotate3d on  % permitir giro con el ratón
% % 4-G  Estadísticas
% for k = 1:3
%     fprintf('Tiempo medio con ≥%d satélites visibles: %.2f h\n', k, mean(vis_hours(:,k)));
% end
% 
% %% Función local: Fibonacci lattice  (¡debe ir al final del script!)
% function pts = fibonacci_sphere(N)
%     phi = (1 + sqrt(5)) / 2;                 % número áureo
%     k   = (0:N-1)';                          % índices 0..N-1
%     y   = 1 - 2*k/(N-1);                     % componente z (-1..1)
%     theta = 2*pi*k*phi;                      % ángulo azimutal
%     r   = sqrt(1 - y.^2);                    % radio en plano XY
%     pts = [r.*cos(theta) , r.*sin(theta) , y]; % (N×3)
% end
% 
