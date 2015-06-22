% Gruppenmitglieder:
% Ben Romdhane, Marouane
% Chihi, Anis
% Ding, Liang
% Gao, Maolin
% Nafouki, Chiraz

%% Bilder laden
close all
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'k',0.03,'N',100,'tile_size',[200,200],'min_dist',23,'segment_length',25,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'k',0.03,'N',100,'tile_size',[200,200],'min_dist',23,'segment_length',25,'do_plot',false);

%% Korrespondenzschätzung
Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2,'min_corr',0.95,'window_length',25,'do_plot',false);

%% Extraktion der Essentiellen Matrix aus robusten Korrespondenzen
% Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus'
%
Korrespondenzen_robust = F_ransac(Korrespondenzen);
size(Korrespondenzen_robust)
%
% Zeigen Sie die robusten Korrespondenzpunkte an
%
showKP(IGray1,IGray2,Korrespondenzen_robust);
%% Berechne die Essentielle Matrix
%
load('K.mat');
E = achtpunktalgorithmus(Korrespondenzen_robust,K);