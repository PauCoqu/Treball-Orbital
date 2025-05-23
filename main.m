filename = 'nom_del_fitxer.txt';
data = readmatrix(filename);


% Obre el fitxer
fid = fopen('data.txt');

% Format: 6 primeres columnes (textuals o mixtes) + 13 columnes numèriques
formatSpec = '%s %d %d %f %s %d %f %f %f %f %f %f %f %f %f %f %f %f %f';

% Llegeix totes les columnes
data = textscan(fid, formatSpec);

% Tanca el fitxer
fclose(fid);