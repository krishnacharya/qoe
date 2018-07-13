function userVec = getUserVec(idx, maxUserVec)
% This function converts the scalar idx to a vector userVec where
% userVec is the state of the system as encoded by idx, i.e., userVec(ii) is 
% the number of users of class ii.

% INPUT
% idx: index which encodes number of users in system
% maxUserVec: Vector (N_1, N_2, ..., N_K) with N_i max. number of class i users allowed in system.

% OUTPUT
% userVec: vector (n_1, n_2, ..., n_K) of users in system (as encoded by idx)

userVec = zeros(1,length(maxUserVec));
idx = idx - 1;
for jj = length(maxUserVec):-1:1,
    remainder = mod(idx,(maxUserVec(jj) + 1));
    quotient = (idx - remainder)/(maxUserVec(jj) + 1);
    idx = quotient;
    userVec(jj) = remainder;
end;
