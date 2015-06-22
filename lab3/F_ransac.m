function [Korrespondenzen_robust]=F_ransac(Korrespondenzen,varargin)
%% inputparser
P           =   inputParser;
P.addRequired('Korrespondenzen');
P.addOptional('epsilon',0.5);
P.addOptional('p',0.99);
P.addOptional('tolerance',0.1);
P.parse(Korrespondenzen,varargin{:});

KPs         =   P.Results.Korrespondenzen;
epsilon     =   P.Results.epsilon;
p           =   P.Results.p;
tao         =   P.Results.tolerance;

%% main part
numKP   =   length(KPs);
k       =   8;
if numKP < k;
    fprintf('Korrespondenzpaaren weniger als 8\n');
    Korrespondenzen_robust = 0;
    return
end

% number of iteration
s       =   log(1-p) / log(1-(1-epsilon)^k);
s       =   ceil(s);

% transform to the homogenous coordinates
x1          =   ones(3,numKP);
x2          =   ones(3,numKP);
x1(1:2,:)   =   KPs(1:2,:);
x2(1:2,:)   =   KPs(3:4,:);

% initial setting
KP_inlier       =   [];     % register of the robust corespondence points
numInlierMax    =   0;      % register for the largest set of consensus
Err             =   0;      % error to be campared for every correspondece pair
e3              =   [0 0 1]';
e3Dach          =   DachVektor(e3);

% random sample consensus
for i = 1 : s
    idx             =   randperm(numKP, k);                 % randomly chosen k points
    F               =   achtpunktalgorithmus(KPs(:, idx));   % computer the foundamental matrix by 8pa
    Sampson_Dis     =   Sampson_Distance(x1, x2, e3Dach, F);    % compute the sampson disntances
    numinlier       =   length( Sampson_Dis( Sampson_Dis < tao)); % compare with the threshold tao
    
    if (numinlier > numInlierMax)   % the new consensus set is larger
                                    % => update the set of consensus
        numInlierMax    =   numinlier;
        KP_inlier       =   KPs(:, Sampson_Dis < tao);
        Err             =   sum( Sampson_Dis ( Sampson_Dis < tao ) );
        
    elseif (numinlier == numInlierMax) % new number of inlier points equals 
                                        % => compare error
        Err_current = sum( Sampson_Dis( Sampson_Dis < tao ) );
        
        if ( Err_current < Err);    % the error of the new consensus set is smaller
                                    % => update the set of consensus
            KP_inlier   =   KPs(:,Sampson_Dis < tao);
            Err         =   Err_current;
        end
    end
                 
end

Korrespondenzen_robust  =   KP_inlier;     % robust KPs

end

%% sub-functions
function A = DachVektor(v)
    % build the skew matrix of a vector v (Dach Operator)
    A      = zeros(3);
    A(1,2) = -v(3);
    A(1,3) = v(2);
    A(2,1) = v(3);
    A(2,3) = -v(1);
    A(3,1) = -v(2);
    A(3,2) = v(1);
end

function d = Sampson_Distance(x1,x2,e3,F)
    % the sampson distance for every KP
    num     =   diag( x2' * F * x1).^2;     
    t1      =   e3 * F * x1;
    den1    =   sum(t1.^2)';
    t2      =   x2' * F * e3;
    den2    =   sum(t2.^2,2);
    d       =   num ./ (den1 + den2);
end
