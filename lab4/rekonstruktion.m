function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

%% initial settings
[~,NrOfCorresp] = size(Korrespondenzen);
% transform pixel coordinates to homogeneous coordinates
zCoord  =   ones(1, NrOfCorresp);
KP1     =   [Korrespondenzen(1:2,:); zCoord];
KP2     =   [Korrespondenzen(3:4,:); zCoord];
% transform pixel coordinates to image coordinates
KP1     =   K \ KP1;
KP2     =   K \ KP2;

%% compute the depth of feature points
% from R1, T1 for camera1
d1c1    =   backtracking3D(KP1, KP2, R1, T1);
% from R1, T1 for camera2
d1c2    =   backtracking3D(KP2, KP1, R1, T1);
% from R2, T2 for camera1
d2c1    =   backtracking3D(KP1, KP2, R2, T2);
% from R2, T2 for camera2
d2c2    =   backtracking3D(KP2, KP1, R2, T2);

%% select the right rotation and transformation
numPos_d1c1     =   positiveDepth(d1c1);
numPos_d1c2     =   positiveDepth(d1c2);
numPos_d2c1     =   positiveDepth(d2c1);
numPos_d2c2     =   positiveDepth(d2c2);

numPos_d1       =   numPos_d1c1 + numPos_d1c2;
numPos_d2       =   numPos_d2c1 + numPos_d2c2;

% compare the number of positive depth, and take the bigger one as R and T
if numPos_d1 > numPos_d2    
    T           =   d(end)*T1; %#ok<COLND>
    R           =   R1;
    lambda1     =   d1c1(1 : end-1);
    lambda2     =   d1c2(1 : end-1);
else
    T           =   d(end)*T2; %#ok<COLND>
    R           =   R2;
    lambda1     =   d2c1(1 : end-1);
    lambda2     =   d2c2(1 : end-1);
end

%% construct the depth lambdas and feature point in space coordinate in Camera1
lambdas     = 	[lambda1; lambda2];
P1          =   [KP1; lambda1];         % !!! multiply each point with corresponding depth

%% plot points in space and camera position


end