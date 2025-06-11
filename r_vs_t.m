function sat_list=r_vs_t(Seconds, SatID, r)
% GRAFICAR_DISTANCIA_GEOCENTRICA Grafica la distancia geocéntrica r frente al tiempo
% para múltiples satélites, diferenciándolos por colores.
%
% Inputs:
%   Seconds - tiempo en segundos desde el inicio del día (vector)
%   SatID   - identificador del satélite (vector)
%   r       - distancia geocéntrica (vector, en metros)

    % Convertir el tiempo a horas
    tiempo_horas = Seconds / 3600;

    % Obtener lista de satélites únicos
    sat_list = unique(SatID);

    % Generar colores distintos
    colores = lines(length(sat_list));

    % Crear figura
    figure;
    hold on; grid on;

    % Graficar cada satélite por separado
    for i = 1:length(sat_list)
        sat_actual = sat_list(i);
        indices = SatID == sat_actual;
        plot(tiempo_horas(indices), r(indices)/1e3, '.', ...
             'Color', colores(i,:), ...
             'DisplayName', ['PRN ' num2str(sat_actual)]);
    end

    % Etiquetas y leyenda
    xlabel('Tiempo (h)');
    ylabel('Distancia geocéntrica r (km)');
    title('Evolución de la distancia geocéntrica de los satélites');
    legend show;
end

