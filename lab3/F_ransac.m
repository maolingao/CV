function [Korrespondenzen_robust]=F_ransac(Korrespondenzen,varargin)
%% inputparser
P = inputParser;
P.addRequired('Korrespondenzen');
P.addOptional('epsilon',0.5);
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
s = log(1-p) / log(1-(1-epsilon)^k);
s = ceil(s);                            % ceiling s

e3 = [0 0 1]';
e3Dach = DachVektor(e3);

% transfer to homogenous coordinate
x1      =   ones(3,Nkp);
x2      =   ones(3,Nkp);
x1(1:2,:) = KPs(1:2,:);
x2(1:2,:) = KPs(3:4,:);

% init setting
numInlierMax   =   0;
KP_inlier =   [];   % index of the robust corespondence points
Err     =   0;
for i = 1 : s
    idx         =   randperm(Nkp, k);              % randomly chosen k points
    F           =   achtpunktalgorithmus(KPs(:,idx)); % call 8-point-algorithm
    Sampson_Dis =   Sampson_Distanz(x1,x2,e3Dach,F);    % compute sampson disntances
    numinlier        =   length(Sampson_Dis(Sampson_Dis < tao)); % compare with threshold
    
    if (numinlier > numInlierMax)   % new number of inlier points is more
        numInlierMax   =   numinlier;
        KP_inlier =   KPs(:,Sampson_Dis < tao);
        Err     =   sum( Sampson_Dis( Sampson_Dis < tao) );
        
    elseif (numinlier == numInlierMax) % new number of inlier points equals 
                                        % => compare error
        Err_current = sum( Sampson_Dis( Sampson_Dis < tao ) );
        
        if ( Err_current < Err); % new error is smaller
            KP_inlier =   KPs(:,Sampson_Dis < tao);
            Err     =   Err_current;
        end
    end
                 
end

Korrespondenzen_robust = KP_inlier;     % robust KPs



end

function A = DachVektor(v)
    % build skew matrix of vector v (Dach Operator)
    A = zeros(3);
    A(1,2) = -v(3);
    A(1,3) = v(2);
    A(2,1) = v(3);
    A(2,3) = -v(1);
    A(3,1) = -v(2);
    A(3,2) = v(1);
end

function d = Sampson_Distanz(x1,x2,e3,F)
    % sampson distance for every KP
    num = diag(x2'* F * x1).^2;     
    t1 = e3 * F * x1;
    den1 = sum(t1.^2)';
    t2 = x2' * F * e3;
    den2 = sum(t2.^2,2);
    d = num ./ (den1 + den2);
end
