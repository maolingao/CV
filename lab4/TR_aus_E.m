function [T1,R1,T2,R2]=TR_aus_E(E)

[U,S,V] = svd(E);

R_zp = [0 -1 0; 1 0 0; 0 0 1];
R_zm = [0 1 0; -1 0 0; 0 0 1];

% Check for U and V if they are rotation matrices
% det(U)
% det(V)
if det(U) < 0
    U = U * [1 0 0; 0 1 0; 0 0 -1];
end

if det(V) < 0
   V = V * [1 0 0; 0 1 0; 0 0 -1];
end
% det(U)
% det(V)
% Solution 1 with R_zp
R1 = U * R_zp' * V';
T1_hat = U * R_zp * S * U';
T1(1,1) = T1_hat(3,2);
T1(2,1) = T1_hat(1,3);
T1(3,1) = T1_hat(2,1);

%T1_hat*R1    % test  == E

% Solution 2 with R_zm
R2 = U * R_zm' * V';
T2_hat = U * R_zm * S * U';
T2(1,1) = T2_hat(3,2);
T2(2,1) = T2_hat(1,3);
T2(3,1) = T2_hat(2,1);

%T2_hat*R2    % test   == E

% det(R1)
% det(R2)
end