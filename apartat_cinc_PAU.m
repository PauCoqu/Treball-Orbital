function [r_ECI, v_ECI, lambda_deg, phi_deg] = apartat_cinc_PAU(SatID, Seconds, DoY, Year, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z)
% APARTAT 5: Transformació del sistema TRF (rotant amb la Terra) al CRF/ECI (inercial)

    %% a) Temps GPS → UTC → UT1
    leap_seconds = 18;
    UT1_UTC = -0.559267;  % IERS Bulletin A, MJD 57103 (22/03/2015)
    UTC_seconds = Seconds - leap_seconds;
    UT1_seconds = UTC_seconds + UT1_UTC;

    %% b) Càlcul del GMST
    omega_E = 7.2921150e-5;   % [rad/s]
    JD0 = 2451545.0;          % J2000.0
    datetime_base = datetime(Year,1,1,0,0,0) + days(DoY - 1);
    JD_UT1 = juliandate(datetime_base) + UT1_seconds / 86400;
    TJC = (JD_UT1 - JD0) / 36525;
    theta_GMST0_sec = 24110.54841 + 8640184.812866*TJC + 0.093104*TJC.^2 - 6.2e-6*TJC.^3;
    theta_GMST_sec = 1.002737909350795 * UT1_seconds + theta_GMST0_sec;
    theta_GMST_sec = mod(theta_GMST_sec, 86400);
    theta_GMST_rad = 2 * pi * theta_GMST_sec / 86400;

    %% c) Rotació TRF → ECI
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

        % d) Correcció velocitat per rotació terrestre
        omega_vec = [0; 0; omega_E];
        v_ECI(:,i) = v_rot + cross(omega_vec, r_rot);
        r_ECI(:,i) = r_rot;
    end

    %% Coordenades esfèriques (λ, φ)
    r_norm = sqrt(sum(r_ECI.^2, 1));
    lambda = atan2(r_ECI(2,:), r_ECI(1,:));
    phi = asin(r_ECI(3,:) ./ r_norm);

    lambda_deg = rad2deg(lambda(:));
    phi_deg = rad2deg(phi(:));

    %% e) Gràfica 3D: òrbites en ECI
    PRNs_dibuixar = [14, 18];
    colors = {'r', 'b', 'g', 'k', 'm', 'c'};

    figure;
    hold on; grid on;
    
    % Dibuixar la Terra com una esfera transparent
    [Xs, Ys, Zs] = sphere(50);
    radio_Terra_km = 6371;
    surf(radio_Terra_km*Xs, radio_Terra_km*Ys, radio_Terra_km*Zs, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.5 0.7 1]);

    % Dibuixar òrbites dels satèl·lits
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

    %% f) Gràfica 2D: projecció longitud - latitud
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
end


