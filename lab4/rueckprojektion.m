function repro_error=rueckprojektion(Korrespondenzen, P1, Image2, T, R, K)

[~,Np1] = size(P1);
KP = Korrespondenzen;

Pi0 =[1 0 0 0;
    0 1 0 0;
    0 0 1 0];

% Berechnen die Projektion der aus Bild 1 rekonstruierten 3D-Punkte in Kamera 2
x2_homo = K * Pi0 * [R, T; zeros(1,3), 1] * [P1;ones(1,Np1)];
x2_x = x2_homo(1,:)./x2_homo(3,:);
x2_y = x2_homo(2,:)./x2_homo(3,:);


currentRepError = zeros(1,Np1);

figure('name','Ruechprojektion')
imshow(Image2)
title(['{\color[rgb]{0 1 0}Robuste Merkmale in Bild 2}, {\color[rgb]{1 0 0}rueckprojizierte Punkte} und der jeweilige ','{\color[rgb]{1 1 1}R?ckprojektionsfehler}'])

hold on
for i=1:Np1
    
    PlotObjKorresp=plot(KP(3,i),KP(4,i),'g*');
    
    PlotObjProj=plot(x2_x,x2_y,'r*');

    x = [x2_x(i), KP(3,i)];
    y = [x2_y(i), KP(4,i)];
    line(x,y,'Color','w');

    xy = [x2_x(i)-KP(3,i);x2_y(i)-KP(4,i)]; 
    currentRepError(1,i) = norm(xy);
    tObj=text((x2_x(i)+KP(3,i))/2,(x2_y(i)+KP(4,i))/2,[' \leftarrow ',num2str(currentRepError(1,i))],'Color','w');

end



repro_error = mean(currentRepError);