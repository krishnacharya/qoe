function badStates = encodeBadStates_balk(vQualMat, throughputVec, maxUsersVec)
% Encode state idx as bad for user of class k if channel rate of class-k user is less than 
% lowest video bit-rate when system state is idx.
% We code what the user sees in the system (state of system is what the user sees)
numClasses = length(maxUsersVec);
numStates = prod(maxUsersVec + 1) + 1;
badStates = zeros(numClasses,numStates);
for idx = 1:numStates,
    userVec = getUserVec(idx, maxUsersVec); % idx --> (i_1, i_2, ..., i_K)
    for kk = 1:numClasses,
        if userVec(kk) < maxUsersVec(kk),
            actualUserVec = userVec; % tagged user does not see itself in system!
            actualUserVec(kk) = userVec(kk) + 1;
            idxPlusOne = codeUserVec(actualUserVec, maxUsersVec);
            if(throughputVec(actualUserVec) < vQualMat(kk,idxPlusOne)) % shouldn't this also check the mimn available bit rate?
                badStates(kk, idx) = 1;
            end
        else
            badStates(kk,idx) = -1;% encode the not possible state with -1
        end;
    end
end;
