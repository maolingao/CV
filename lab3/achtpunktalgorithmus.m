function [EF]=achtpunktalgorithmus(Korrespondenzen,K)

%% general rearrange coordinates
% homogeneous coordinates in pixel
[~,numKp]   =   size(Korrespondenzen);
z           =   ones(1,numKp);
kp1         =   [Korrespondenzen(1:2,:);z];
kp2         =   [Korrespondenzen(3:4,:);z];

% homogeneous coordinates in pixel to homogeneous coordinates in image
if exist('K','var')
    % if intrinsic parameter of camera is available
    % => be able to calculate E matrix.
    kp1     =   K \ kp1;                                % calibrated coordi.
    kp2     =   K \ kp2;
end

%% build matrix A
A_reg       =   kron(kp1,kp2);
A           =   A_reg(:, 1 : numKp+1 : end);
A           =   A';

%% SVD of A
[~,~,V_a]   =   svd(A);
G_s         =   V_a(:,9);                               % last column of rotation matrix V
G           =   reshape(G_s,3,3);

%% SVD of G
[U_g,S_g,V_g] = svd(G);

%% output corresponding matrix E/F
if exist('K','var')
    sigma   =   1;
    EF      =   U_g * diag([sigma,sigma,0]) * V_g';     % essential matrix E
else
    S_g(3,3)    =   0;
    EF          =   U_g * S_g * V_g';                   % fundamental matrix F
end
end
