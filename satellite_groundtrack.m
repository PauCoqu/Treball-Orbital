function satellite_groundtrack(lambda_deg, phi_deg, SatID)
% GRAFICAR_GROUND_TRACKS_PUNTOS - Muestra los ground tracks de GSAT201 y GSAT202
% como puntos sobre un mapa de fondo (sin toolboxes especiales).
%
% Inputs:
%   lambda_deg - vector de longitudes (en grados)
%   phi_deg    - vector de latitudes (en grados)
%   SatID      - vector con el identificador del satélite (PRN)
%
% Usa una imagen de mapa de fondo 'World_map.png' en la carpeta actual.

    % Crear figura
    figure;

    % Cargar imagen del mapa del mundo (ajusta el nombre si es necesario)
    img = imread('World_map.png'); % Mapa en proyección equirectangular
    image([-180 180], [-90 90], flipud(img)); % Ajustar ejes
    set(gca,'YDir','normal'); % Orientación correcta del eje Y
    hold on;

    % Ajustar límites de los ejes
    xlim([-180 180]);
    ylim([-90 90]);

    % Dibujar ground tracks como puntos
    plot(lambda_deg(SatID == 14), phi_deg(SatID == 14), 'r.', ...
         'LineWidth', 1.5, 'DisplayName', 'GSAT201 (PRN 14)');
    plot(lambda_deg(SatID == 18), phi_deg(SatID == 18), 'b.', ...
         'LineWidth', 1.5, 'DisplayName', 'GSAT202 (PRN 18)');

    % Etiquetas y título
    xlabel('Longitud (°)');
    ylabel('Latitud (°)');
    title('Ground tracks de GSAT201 y GSAT202 (TRF)');
    legend show;
end
