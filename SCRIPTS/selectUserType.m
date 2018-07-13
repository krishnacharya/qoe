function idx = selectUserType(rateVec)
% Given a vector of weights, generate idx randomly with probability (w.p.) proportional to the weights.
% Suppose rateVec = [0.2, 0.6, 0.2] then output idx = 1 w.p. 0.2, 2 w.p. 0.6 and 3 w.p. 0.2
idVec = 1:length(rateVec);
if (sum(rateVec) > 0)
    idx = min( idVec( rand() < cumsum(rateVec)/sum(rateVec) ) );
else
    idx = 0
end
end
