function v_vs_t(SatID, x, y, z, v_x, v_y, v_z)
% GRAFICAR_VELOCIDAD_VS_POSICION Graficar velocidad total vs. distancia geocéntrica
%
% Inputs:
%   SatID - ID del satélite
%   x, y, z - coordenadas de posición
%   v_x, v_y, v_z - componentes de velocidad

    % Calcular magnitudes
    r = sqrt(x.^2 + y.^2 + z.^2);
    v = sqrt(v_x.^2 + v_y.^2 + v_z.^2);

    % Lista de satélites únicos
    sat_list = unique(SatID);
    colores = lines(length(sat_list));

    figure; hold on; grid on;
    for i = 1:length(sat_list)
        indices = SatID == sat_list(i);
        plot(r(indices)/1e3, v(indices)/1e3, '.', ...
             'Color', colores(i,:), ...
             'DisplayName', ['PRN ' num2str(sat_list(i))]);
    end

    xlabel('Distancia geocéntrica r (km)');
    ylabel('Velocidad orbital v (km/s)');
    title('Velocidad vs. Posición (r) de los satélites');
    legend show;
end

