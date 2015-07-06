function [T1,R1,T2,R2]=TR_aus_E(E)
% In dieser Funktion sollen die moegleichen euklidischen Transformationen
% aus der essentiellen Matrix extrahiert werden
RzP = [0 -1 0; 1 0 0; 0 0 1]; % mit positiv pi Haefte
RzN = [0 1 0; -1 0 0; 0 0 1]; % mit negativ pi Haefte
[U Sigma V] = svd(E); % Singulaerwertzerlegung der essentiellen Matrix

% Loesung 1 mit RzP
R1 = U * RzP' * V';
T1_Dach = U * RzP * Sigma * U';
T1 = SkewMatrix2Vector(T1_Dach);

% Loesung 2 mit RzN
R2 = U * RzN' * V';
T2_Dach = U * RzN * Sigma * U';
T2 = SkewMatrix2Vector(T2_Dach);

end

function v = SkewMatrix2Vector(A)

v(1) = A(3,2);
v(2) = A(1,3);
v(3) = A(2,1);

v = v';
end