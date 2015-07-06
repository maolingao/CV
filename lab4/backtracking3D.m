function d = backtracking3D(KP0,KP1,R,T,aux)
% compute the depth information of the 3d points and the scale of the transformation.
% input : correspondence points KP0 (the reference of the depth
% information) and KP1, rotation R, transformation T.
% output : vector d = [lamdas; gamma];

% construct matrix M
M   =   buildMatrixM(KP0,KP1,R,T,aux);

% solve linear system M*d = 0
d   =   linearMinimizer(M,aux);
end

function [M] = buildMatrixM(KP0,KP1,R,T,aux)
% build matrix M for linear system M*d=0

%# some intermediate registers
KP0_cell    =   num2cell(KP0, 1);
RKP0        =   cellfun(@(a)R*a, KP0_cell, 'UniformOutput', 0); % rotated KP0

%# apply hat-operator on KP1
KP1_cell    =   num2cell(KP1, 1);
KP1_hat     =   cellfun(@skew_matrix, KP1_cell, 'UniformOutput', 0);

%# build M
Rpart       =   cellfun(@(a,b)a*b, KP1_hat, RKP0, 'UniformOutput', 0);
Tpart       =   cellfun(@(a)a*T, KP1_hat, 'UniformOutput', 0);

MTpart      =   cell2mat(Tpart');                       % last column of M
MRpart      =   zeros(length(MTpart),length(Tpart));
for k = 1 : length(Tpart)                               % rotation part of M 
    MRpart(1 + 3*(k-1) : 3*k, k)  =   Rpart{k};
end

M   =   [MRpart MTpart];                                % matrix M

end

function [d] = linearMinimizer(M,aux) %#ok<*INUSD>
% solve linear system M*d=0

[~,~,V]     =   svd(M);             % svd of M
d           =   V(:,end);           % the last column of V as the solution
end

function [Vhat]=skew_matrix(V)
% hat-operator
Vhat        =   [0 -V(3) V(2); V(3) 0 -V(1); -V(2) V(1) 0];
end