function [r, lambda_deg, phi_deg, lambda, phi] = Cartesianes_geocentriques(x, y, z)

    % Módulo del vector posición
    r = sqrt(x.^2 + y.^2 + z.^2);

    % Longitud (lambda), con atan2 para tener el ángulo en el cuadrante correcto
    lambda = atan2(y, x);           % En radianes

    % Latitud geocéntrica (phi)
    phi = asin(z ./ r);            % En radianes

    % Convertir a grados
    lambda_deg = rad2deg(lambda);
    phi_deg = rad2deg(phi);
end
