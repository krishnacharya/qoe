function pi_W = getStationaryDist(TPM)
% Given the MC Transition probability matrix TPM, evaluate the stattionary
% probabilities
%Usage: pi_W = getStationaryDist(TPM)

if(length(TPM) < 5000),
    opts.disp = 0;
    [V, D] = eigs(double(TPM'),1, 'LM', opts); % finds left-evec of TPM corr. to eval = 1
else % for large Markov chains, try iterating pi(k+1) = pi(k)*P
    [V, D] = iterator(TPM);
end;
if(abs(D - 1) > 1e-5),
    warning(0, 'GetStationaryDist_W_D.m: Check Transition Matrix');
    exit;
end;
V = abs(V);
pi_W = V/sum(V); % Normalize to get probability distribution