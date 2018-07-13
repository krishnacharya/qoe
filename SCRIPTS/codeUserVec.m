function idx = codeUserVec(userVec, maxUsersVec)
% Convert vector of users to a scalar index
% f: (n_1,n_2, ..., n_K) --> idx
% INPUT:
% userVec: vector of users (n_1, n_2, ..., n_K)
% maxUsersVec: Vector of max. number of users permitted in system for each class
% OUTPUT:
% idx: encodes vector of users into a scalar
idx = userVec(length(maxUsersVec)) + 1;
prod = 1;
for jj=length(maxUsersVec):-1:2
    prod = prod*(maxUsersVec(jj)+1);
    idx = idx + userVec(jj-1)*prod;
end
