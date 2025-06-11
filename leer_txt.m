function [Year, DoY, Seconds, Constellation, SatID, x_TRF, y_TRF, z_TRF, v_x, v_y, v_z, clock_offset] = leer_txt (filename);
% Obre el fitxer
fid = fopen(filename, 'r');
% Format del fitxer
formatSpec = '%s %f %f %f %s %f %f %f %f %f %f %f %f'; %s = text ; d= dades sense comes ; f = dades amb comes
% Llegeix totes les columnes
data = textscan(fid, formatSpec);

% Tanca el fitxer
fclose(fid);

% Assigna cada columna a una variable separada (textual o numèrica)
%col1_tag = data{1};   % 'SATPVT' --> ens la pela
Year = data{2};   % 2015
DoY = data{3};   % 081  (DoY = Day of Year)
Seconds = data{4};   % 0.00, 300.00... (Seconds elapsed since the beginning of the day in GPS time.)
Constellation = data{5};   % 'GAL'
SatID = data{6};   % 11, 12, 14... 

% Les 13 columnes de dades numèriques
x_TRF = data{7};   % columna 7 (x_TRF)
y_TRF = data{8};   % columna 8
z_TRF = data{9};   % columna 9
v_x = data{10};  % columna 10
v_y = data{11};  % columna 11
v_z = data{12};  % columna 12
clock_offset = data{13};  % columna 13

end