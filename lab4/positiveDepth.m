function numPosDepth = positiveDepth(d, aux)
% count the number of positive depth in d = [lamdas; gamma];

lambda      = d(1:end-1);
numPosDepth = sum(lambda>0);

end