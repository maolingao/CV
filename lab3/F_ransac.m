function [Korrespondenzen_robust]=F_ransac(Korrespondenzen,varargin)
%% inputparser
P = inputParser;
P.addRequired('Korrespondenzen');
P.addOptional('epsilon',0.2);
P.addOptional('p',0.99);
P.addOptional('tolerance',0.1);
P.parse(Korrespondenzen,varargin{:});


KPs = P.Results.Korrespondenzen;
epsilon = P.Results.epsilon;
p = P.Results.p;
tao= P.Results.tolerance;
%%

Nkp = length(KPs);
k = 8;
if Nkp < k;
    fprintf('Korrespondenzpaaren weniger als 8\n');
    Korrespondenzen_robust = 0;
    return
end

% Iterationszahl berechnen
s = log(1-p)/log(1-(1-epsilon)^k);

e3 = [0 0 1]';
e3Dach = DachVektor(e3);

x1 = ones(3,Nkp);
x2 = ones(3,Nkp);
x1(1:2,:) = KPs(1:2,:);
x2(1:2,:) = KPs(3:4,:);

N_in = 0;
N_max = 0;
KP_need = [];
Err = 0;
for i=1:s
    % zufaellige k Punkte
    idx = randperm(Nkp,k);
    F = achtpunktalgorithmus(KPs(:,idx));
    Sampson_Dis = Sampson_Distanz(x1,x2,e3Dach,F);
    N_in = length(Sampson_Dis(Sampson_Dis < tao));
    if(N_in > N_max)
        N_max = N_in;
        KP_need = KPs(:,Sampson_Dis<tao);
        Err = sum(Sampson_Dis(Sampson_Dis<tao));
    elseif(N_in == N_max)
        if(Err>sum(Sampson_Dis(Sampson_Dis<tao)));
            N_max = N_in;
            KP_need = KPs(:,Sampson_Dis<tao);
            Err = sum(Sampson_Dis(Sampson_Dis<tao));
        end
    end
                 
end
keyboard
Korrespondenzen_robust = KPs(:,KP_need);



end

function A = DachVektor(v)
A = zeros(3);
A(1,2) = -v(3);
A(1,3) = v(2);
A(2,1) = v(3);
A(2,3) = -v(1);
A(3,1) = -v(2);
A(3,2) = v(1);
end

function d = Sampson_Distanz(x1,x2,e3,F)

num = diag(x2'* F * x1).^2;
t1 = e3*F*x1;
den1 = diag(t1'*t1);
t2 = x2'*F*e3;
den2 = diag(t2*t2');
d = num./(den1+den2);

end
