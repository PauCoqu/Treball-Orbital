function [x_CRF, y_CRF, z_CRF, v_x_CRF, v_y_CRF, v_z_CRF] = trf_to_crf(Year_UT, DoY_UT, Seconds_UT, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z)
%TRF_TO_CRF Convert Earth fixed coordinates to Celestial Reference Frame.
%   [x_CRF, y_CRF, z_CRF, v_x_CRF, v_y_CRF, v_z_CRF] = TRF_TO_CRF(Year_UT,
%   DoY_UT, Seconds_UT, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z) rotates the
%   position and velocity vectors from the terrestrial reference frame (TRF)
%   to the celestial reference frame (CRF) using the UT1 time. The UT1 time
%   is given by year, day-of-year and seconds of day.
%
%   The Greenwich mean sidereal time (GMST) is computed with the IAU 1982
%   conventions and the vectors are rotated about the Z axis by -GMST.

    N = numel(Seconds_UT);
    x_CRF   = zeros(size(x_TRF));
    y_CRF   = zeros(size(y_TRF));
    z_CRF   = zeros(size(z_TRF));
    v_x_CRF = zeros(size(v_x));
    v_y_CRF = zeros(size(v_y));
    v_z_CRF = zeros(size(v_z));

    for k = 1:N
        % Convert UT1 to Julian Date
        dt = datetime(Year_UT(k),1,1) + days(DoY_UT(k)-1) + seconds(Seconds_UT(k));
        JD = juliandate(dt);

        % Julian centuries from J2000
        T = (JD - 2451545.0) / 36525;

        % GMST from IAU 1982 (in seconds)
        theta_sec = 67310.54841 + (876600*3600 + 8640184.812866)*T ...
                     + 0.093104*T^2 - 6.2e-6*T^3;
        % Convert to radians in [0,2*pi)
        theta = mod(theta_sec, 86400) * (2*pi/86400);

        % Rotation matrix about Z axis
        R = [ cos(theta)  sin(theta)  0;\
              -sin(theta)  cos(theta)  0;\
                    0          0      1];

        % Rotate position and velocity
        r_crf = R * [x_TRF(k); y_TRF(k); z_TRF(k)];
        v_crf = R * [v_x(k);  v_y(k);  v_z(k)];

        x_CRF(k)   = r_crf(1);
        y_CRF(k)   = r_crf(2);
        z_CRF(k)   = r_crf(3);
        v_x_CRF(k) = v_crf(1);
        v_y_CRF(k) = v_crf(2);
        v_z_CRF(k) = v_crf(3);
    end
end
