function ENU = xyz_2_ENU(x_obs, y_obs, z_obs, lambda_rad4, phi_rad4, N)

    for i = 1:N
        R_matriu = [
            -sin(lambda_rad4(i)),            cos(lambda_rad4(i)),            0;
            -cos(lambda_rad4(i))*sin(phi_rad4(i)),  -sin(lambda_rad4(i))*sin(phi_rad4(i)),   cos(phi_rad4(i));
             cos(lambda_rad4(i))*cos(phi_rad4(i)),   sin(lambda_rad4(i))*cos(phi_rad4(i)),   sin(phi_rad4(i))
        ];
    
    ENU(i, :) = (R_matriu * [x_obs(i); y_obs(i); z_obs(i)])';

    end
end

