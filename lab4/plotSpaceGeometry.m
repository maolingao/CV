function plotSpaceGeometry(PointInSpace1,R,T,aux)
% show 3d point position in space and camera position before and after
% rigid body motion

figure(100)
title('3D-Punkte bezüglich des Kamerasystems 1')
xlabel('x-Achse')
ylabel('y-Achse')
zlabel('z-Achse')

xmin = min(PointInSpace1(1,:))-0.3;
xmax = max(PointInSpace1(1,:))+0.3;
ymin = min(PointInSpace1(2,:))-0.3;
ymax = max(PointInSpace1(2,:))+0.3;
zmin = 0;
zmax = max(PointInSpace1(3,:))+0.5;

axis([xmin xmax ymin ymax zmin zmax]) 

hold on

for k=1 : size(PointInSpace1,2)   
    currentColor = rand(1,3);
    plot3(PointInSpace1(1,k),PointInSpace1(1,k),PointInSpace1(3,k),'marker','o','Color',currentColor)
    text(PointInSpace1(1,k),PointInSpace1(1,k),PointInSpace1(3,k),num2str(k))
    %plot3(PointInSpace2(1,k),PointInSpace2(2,k),PointInSpace2(3,k),'marker','+','Color',currentColor)
    %text(PointInSpace2(1,k),PointInSpace2(2,k),PointInSpace2(3,k),num2str(k))
    grid on
end

cam1 = [0;0;1];
cam2 = R*cam1+T;

plot3(cam1(1),cam1(2),cam1(3),'marker','s','Color','k')
text(cam1(1),cam1(2),cam1(3),'Kamera 1')
plot3(cam2(1),cam2(2),cam2(3),'marker','s','Color','r')
text(cam2(1),cam2(2),cam2(3),'Kamera 2')

az = 20;
el = 60;
view(az, el);
end