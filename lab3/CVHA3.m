% Gruppenmitglieder:

%% Bilder laden
% Image1 = imread('szeneL.jpg');
% IGray1 = rgb_to_gray(Image1);
% 
% Image2 = imread('szeneR.jpg');
% IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
% Merkmale1 = harris_detektor(IGray1,'k',0.03,'N',100,'min_dist',20,'segment_length',9,'do_plot',false);
% Merkmale2 = harris_detektor(IGray2,'k',0.03,'N',100,'min_dist',20,'segment_length',9,'do_plot',false);

%% Korrespondenzschätzung
% Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2,'do_plot',false);
% save('Korrespondenzen.mat','Korrespondenzen')
load('Korrespondenzen.mat')
%% Extraktion der Essentiellen Matrix aus robusten Korrespondenzen
% Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus'
%
Korrespondenzen_robust = F_ransac(Korrespondenzen);
%
% Zeigen Sie die robusten Korrespondenzpunkte an

%% Berechne die Essentielle Matrix
%
F = achtpunktalgorithmus(Korrespondenzen_robust);
load('K.mat');
E = achtpunktalgorithmus(Korrespondenzen_robust,K);
inv(K')*E*inv(K) - F