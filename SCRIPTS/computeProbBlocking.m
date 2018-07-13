function probBlocking = computeProbBlocking(pi, maxUsersVec)
% return probability of users seeing N_i users in the system i.e., prob. of user dropping per class
% for e.g for class 1 the blocking states are (N1,0,0...), (N1,0,1,etc..)
numClasses = length(maxUsersVec);
numStates = prod(maxUsersVec + 1);
probBlocking = zeros(1,numClasses);
for jj = 1:numStates,
    userVec = getUserVec(jj, maxUsersVec); % idx --> (i_1, i_2, ..., i_K)
    for kk = 1:numClasses,
        if(userVec(kk) == maxUsersVec(kk))
            probBlocking(kk) = probBlocking(kk) + pi(jj);
        end
    end;
end;
