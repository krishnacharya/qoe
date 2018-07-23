function pi = getStationaryDistDirect(TRM)
% Given the MC Transition probability matrix TPM, evaluate the stattionary
% probabilities
%Usage: pi_W = getStationaryDist(TPM), Note we use matrix inversion complexity is O(N^3)
N = length(TRM);
pi = ones(1, N) * inv(TRM + ones(N,N));
