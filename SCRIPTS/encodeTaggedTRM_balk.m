function trm = encodeTaggedTRM_balk(arrivalRateVec, averageVideoSizeVec, ...
    gammaVec, thresVec, vQualMat, throughputVec, maxUsersVec, classId)
% The code gives the TRM corresponding to tagged user Markov chain
% Here we code what the user sees in the system (state of system is what the user sees)
% 
numClasses = length(maxUsersVec);
numStates = prod(maxUsersVec + 1); 
trm = zeros(numStates + 1);% additional state for user balking
for idx = 1:numStates,
    userVec = getUserVec(idx, maxUsersVec); % idx --> (i_1, i_2, ..., i_K)
    if(userVec(classId) == maxUsersVec(classId))
        continue;
    end;
    
    % Encode arrival of each class
    for ii=1:numClasses,
        nextUserVec = userVec;
        if(ii ~= classId) % for other user classes
            if userVec(ii) < maxUsersVec(ii),
                nextUserVec(ii) = userVec(ii) + 1;
                nextIdx = codeUserVec(nextUserVec, maxUsersVec);
                trm(idx, nextIdx) = arrivalRateVec(ii);
            end
        else % for tagged user class ***************IMP************
            if userVec(ii) < maxUsersVec(ii) - 1,
                nextUserVec(ii) = userVec(ii) + 1;
                nextIdx = codeUserVec(nextUserVec, maxUsersVec);
                trm(idx, nextIdx) = arrivalRateVec(ii);
            end            
        end
    end
    
    % Encode departure of each class
    for ii=1:numClasses,
        nextUserVec = userVec;
        if(ii ~= classId) % for other user classes
            if userVec(ii) > 0,
                nextUserVec(ii) = userVec(ii) - 1;
                nextIdx = codeUserVec(nextUserVec, maxUsersVec);
                trm(idx, nextIdx) = trm(idx, nextIdx) + userVec(ii) * throughputVec(userVec(ii)) / (vQualMat(ii,idx) * averageVideoSizeVec(ii));
            end
        else % for tagged user class
            if userVec(ii) < maxUsersVec(ii) && userVec(ii) > 0,
                nextUserVec(ii) = userVec(ii) - 1;
                nextIdx = codeUserVec(nextUserVec, maxUsersVec);
                
                actualUserVec = userVec; % tagged user does not see itself in system!
                actualUserVec(ii) = userVec(ii) + 1;
                idxPlusOne = codeUserVec(actualUserVec, maxUsersVec);        
                trm(idx, nextIdx) = trm(idx, nextIdx) + userVec(ii) * throughputVec(actualUserVec) / (vQualMat(ii, idxPlusOne) * averageVideoSizeVec(ii));
            end
        end
    end
    
    % Encode balking of each class
    for ii=1:numClasses,
        nextUserVec = userVec;
        if(ii ~= classId) % for other user classes
            if userVec(ii) > 0,
                nextUserVec(ii) = userVec(ii) - 1;
                nextIdx = codeUserVec(nextUserVec, maxUsersVec);
                trm(idx, nextIdx) = trm(idx, nextIdx) + userVec(ii)* ...
                    (vQualMat(ii,idx) < thresVec(ii)) ...
                    * (gammaVec(ii) /vQualMat(ii,idx));
            end
        else % for tagged user class
            if userVec(ii) < maxUsersVec(ii) && userVec(ii) > 0,
                nextUserVec(ii) = userVec(ii) - 1;
                nextIdx = codeUserVec(nextUserVec, maxUsersVec);
                
                actualUserVec = userVec; % tagged user does not see itself in system!
                actualUserVec(ii) = userVec(ii) + 1;
                idxPlusOne = codeUserVec(actualUserVec, maxUsersVec);
                trm(idx, nextIdx) = trm(idx, nextIdx) + userVec(ii)* ... % tagged user does not see itself in system!
                    (vQualMat(ii,idxPlusOne) < thresVec(ii)) ...
                    * (gammaVec(ii) / vQualMat(ii,idxPlusOne));
            end;
        end
    end
    
    
    % Encode tagged user departure ( absorbing state (N_1, N_2, ...; N_K) )
    if userVec(classId) < maxUsersVec(classId),
        nextUserVec = maxUsersVec;
        nextIdx = codeUserVec(nextUserVec, maxUsersVec);
        actualUserVec = userVec; % tagged user does not see itself in system!
        actualUserVec(classId) = userVec(classId) + 1;
        idxPlusOne = codeUserVec(actualUserVec, maxUsersVec);
        trm(idx, nextIdx) = trm(idx, nextIdx) + ... % tagged user does not see itself in system!
            (throughputVec(actualUserVec) / vQualMat(classId,idxPlusOne)) / averageVideoSizeVec(classId);
    end;
    
    % Encode tagged user balking (absorbing state (N_1, N_2, ...; N_K) + 1)
    if userVec(classId) < maxUsersVec(classId),
        nextUserVec = maxUsersVec;
        nextIdx = codeUserVec(nextUserVec, maxUsersVec) + 1;
        actualUserVec = userVec; % tagged user does not see itself in system!
        actualUserVec(classId) = userVec(classId) + 1;
        idxPlusOne = codeUserVec(actualUserVec, maxUsersVec);
        trm(idx, nextIdx) = trm(idx, nextIdx) + ... % tagged user does not see itself in system!
            (vQualMat(classId,idxPlusOne) < thresVec(classId)) ...
            * (gammaVec(classId) /vQualMat(classId,idxPlusOne));
    end;
    
    
    trm(idx,idx) = -sum(trm(idx,:));
end
