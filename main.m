%% CODI PROYECTO ORBITAL
%---------------------------------------------
% Fecha: 23/05/2025
% Miembros: Antonio Luque, Pau Cornudella i Alex Aleñà
% Asignatura: AERODINÁMICA, MECÁNICA DE VUELO Y ORBITAL
% Grupo: 6
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

%% 5. Coordinate transformation from TRF to ECI
    % a) Temps GPS -> UTC -> UT1
    leap_seconds = 18;
    UT1_UTC = -0.559267;  % IERS Bulletin A, MJD 57103 (22/03/2015)
    UTC_seconds = Seconds - leap_seconds;
    UT1_seconds = UTC_seconds + UT1_UTC;

    % b) Càlcul del GMST
    omega_E = 7.2921150e-5;   % [rad/s]
    JD0 = 2451545.0;          % J2000.0
    datetime_base = datetime(Year,1,1,0,0,0) + days(DoY - 1);
    JD_UT1 = juliandate(datetime_base) + UT1_seconds / 86400;
    TJC = (JD_UT1 - JD0) / 36525;
    theta_GMST0_sec = 24110.54841 + 8640184.812866*TJC + 0.093104*TJC.^2 - 6.2e-6*TJC.^3;
    theta_GMST_sec = 1.002737909350795 * UT1_seconds + theta_GMST0_sec;
    theta_GMST_sec = mod(theta_GMST_sec, 86400);
    theta_GMST_rad = 2 * pi * theta_GMST_sec / 86400;

    % c) Rotación TRF -> ECI
    R3 = @(theta) [cos(theta) sin(theta) 0;
                  -sin(theta) cos(theta) 0;
                   0          0          1];

    N = length(Seconds);
    r_ECI = zeros(3,N);
    v_ECI = zeros(3,N);

    for i = 1:N
        R = R3(-theta_GMST_rad(i));
        r_TRF = [x_TRF(i); y_TRF(i); z_TRF(i)];
        v_TRF = [v_x(i); v_y(i); v_z(i)];

        r_rot = R * r_TRF;
        v_rot = R * v_TRF;

        % d) Correccion velocidad per rotacion terrestre
        omega_vec = [0; 0; omega_E];
        v_ECI(:,i) = v_rot + cross(omega_vec, r_rot);
        r_ECI(:,i) = r_rot;
    end

    % Coordenadas esfericas (lambda, phi)
    r_norm = sqrt(sum(r_ECI.^2, 1));
    lambda = atan2(r_ECI(2,:), r_ECI(1,:));
    phi = asin(r_ECI(3,:) ./ r_norm);

    lambda_deg = rad2deg(lambda(:));
    phi_deg = rad2deg(phi(:));

    % e) Graficas
    PRNs_dibuixar = [14, 18];
    colors = {'r', 'b', 'g', 'k', 'm', 'c'};

    figure;
    hold on; grid on;
    
    % Esfera Tierra
    [Xs, Ys, Zs] = sphere(50);
    radio_Terra_km = 6371;
    surf(radio_Terra_km*Xs, radio_Terra_km*Ys, radio_Terra_km*Zs, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.5 0.7 1]);

    for k = 1:length(PRNs_dibuixar)
        prn = PRNs_dibuixar(k);
        idx = (SatID == prn);
        if any(idx)
            plot3(r_ECI(1,idx)/1e3, r_ECI(2,idx)/1e3, r_ECI(3,idx)/1e3, ...
                [colors{k} '-'], 'LineWidth', 1.2, 'DisplayName', ['PRN ' num2str(prn)]);
        end
    end
    xlabel('X_{ECI} [km]');
    ylabel('Y_{ECI} [km]');
    zlabel('Z_{ECI} [km]');
    title('Orbitas de los satelites Galileo en sistema ECI');
    legend show;
    axis equal;

    % Gráfica 2D
    figure;
    hold on; grid on;
    for k = 1:length(PRNs_dibuixar)
        prn = PRNs_dibuixar(k);
        idx = (SatID == prn);
        if any(idx)
            plot(lambda_deg(idx), phi_deg(idx), ...
                '.', 'Color', colors{k}, 'DisplayName', ['PRN ' num2str(prn)]);
        end
    end
    xlabel('Longitud (°)');
    ylabel('Latitud (°)');
    title('Ground tracks (2D) en sistema ECI');
    legend show;
    xlim([-180 180]);
    ylim([-90 90]);


%% 6. Osculating orbital elements y gráfica de órbitas en ECI

mu = 3.986004418e14;  % [m^3/s^2]
prns = unique(SatID,'stable');
Nsat = numel(prns);
a_km     = nan(Nsat,1);
e        = nan(Nsat,1);
i_deg    = nan(Nsat,1);
Omega_deg= nan(Nsat,1);
omega_deg= nan(Nsat,1);
nu_deg   = nan(Nsat,1);

for k = 1:Nsat
    % Fórmulas p100,101
    prn = prns(k);
    idx = find(SatID==prn,1,'first');
    % r, v en ECI
    r = r_ECI(:,idx);
    v = v_ECI(:,idx);

    % 1) Distancia r
    R   = norm(r);          % r = √(X^2+Y^2+Z^2)

    % 2) Velocidad escalar v
    V2  = dot(v,v);         

    % 5) Momento angular h
    h     = cross(r,v);
    h_norm= norm(h);

    % 6) Semieje mayor a
    energy = V2/2 - mu/R;   % E = v^2/2 – μ/r
    a      = -mu/(2*energy);
    a_km(k)= a/1e3;

    % 4) Excentricidad e 
    e_vec = cross(v,h)/mu - r/R;
    e(k)  = norm(e_vec);

    % 7) Inclinación i
    i_deg(k) = acosd( h(3)/h_norm );

    % 8) Línea de nodos N
    Nvec  = cross([0;0;1], h);
    Nnorm = norm(Nvec);

    % 9) Ascensión recta Omega
    Om = acosd( Nvec(1)/Nnorm );
    if Nvec(2) < 0, Om = 360-Om; end
    Omega_deg(k) = Om;

    % 10) Argumento del periapsis omega
    om = acosd( dot(Nvec,e_vec)/(Nnorm*e(k)));
    if e_vec(3) < 0, om = 360-om; end
    omega_deg(k) = om;

    % 11) Anomalía verdadera ν 
    nu = acosd( dot(e_vec,r)/(e(k)*R) );
    if dot(r,v) < 0, nu = 360-nu; end
    nu_deg(k) = nu;

end


% Mostrar tabla en la consola
T = table(prns, a_km, e, i_deg, Omega_deg, omega_deg, nu_deg, ...
    'VariableNames',{'PRN','a_km','e','i_deg','Omega_deg','omega_deg','nu_deg'});
disp(T)

figure('Color','w');
hold on; grid on; axis equal tight
view(40,25)

% Esfera de la Tierra
[xe,ye,ze] = sphere(100);
radio_Terra_km = 6371;
    surf(radio_Terra_km*Xs, radio_Terra_km*Ys, radio_Terra_km*Zs, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.5 0.7 1]);

prns = unique(SatID,'stable');
N    = numel(prns);
cols = lines(N);                

% Orbita en km
 for k = 1:N
     idx = SatID==prns(k);
     X = r_ECI(1,idx)/1e3;
     Y = r_ECI(2,idx)/1e3;
    Z = r_ECI(3,idx)/1e3;
     plot3(X, Y, Z, 'LineWidth',1.2, 'Color', cols(k,:), ...
           'DisplayName', ['PRN ' num2str(prns(k))]);
end

% prns_plot = [14, 18];  % PRNs a mostrar
% cols = lines(numel(prns_plot));
% for k = 1:numel(prns_plot)  %14 i 18 només (apartado 6b)
%     idx = SatID == prns_plot(k);
%     X = r_ECI(1,idx)/1e3;
%     Y = r_ECI(2,idx)/1e3;
%     Z = r_ECI(3,idx)/1e3;
%     plot3(X, Y, Z, 'LineWidth',1.2, 'Color', cols(k,:), ...
%           'DisplayName', ['PRN ' num2str(prns_plot(k))]);
%end

xlabel('X_{ECI} (km)')
ylabel('Y_{ECI} (km)')
zlabel('Z_{ECI} (km)')
title('Órbitas de los satélites Galileo en sistema ECI')
legend('Location','northeastoutside')
