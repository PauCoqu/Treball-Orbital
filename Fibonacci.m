function pts = Fibonacci(N)


    phi   = (1 + sqrt(5)) / 2;        % número áureo
    k     = (0:N-1)';                 % índices 0..N-1
    z     = 1 - 2*k/(N-1);            % componente z en [-1,1]
    theta = 2*pi*k/phi;               % ángulo azimutal
    r     = sqrt(1 - z.^2);           % radio en el plano XY

    pts   = [r .* cos(theta), ...     % coordenada x
             r .* sin(theta), ...     % coordenada y
             z];                      % coordenada z
end
