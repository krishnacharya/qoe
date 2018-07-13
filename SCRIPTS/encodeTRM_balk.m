function trm = encodeTRM_balk(arrivalRateVec, averageVideoSizeVec, gammaVec, ...
    thresVec, vQualMat, rateMat, maxUsersVec)
% Encode the trm for the system
% INPUT:
% rateMat: Matrix of average channel rates as seen by users in different states
% 	rateMat(ii,idx) : channel rate of user of class ii when state is idx
% vQualMat: Matrix of average video bit-rate
%	vQualMat(ii,idx): Avg. video bit-rate of user of class ii when state is idx
% IMP: rateMat and vQualMat are identical when l_{min} < r < l_{max}

numStates = prod(maxUsersVec + 1);
trm = zeros(numStates);
for idx = 1:numStates,
    userVec = getUserVec(idx, maxUsersVec); % idx --> (i_1, i_2, ..., i_K)
    % Encode arrival of each class
    for ii=1:length(maxUsersVec),
        nextUserVec = userVec;
        if userVec(ii) < maxUsersVec(ii),
            nextUserVec(ii) = userVec(ii) + 1;
            nextIdx = codeUserVec(nextUserVec, maxUsersVec);
            trm(idx, nextIdx) = arrivalRateVec(ii); %% same arrival rate for a given class
        end
    end
    % Encode departure and user balking for each class
    for ii=1:length(maxUsersVec),
        nextUserVec = userVec;
        if userVec(ii) > 0,
            nextUserVec(ii) = userVec(ii) - 1;
            nextIdx = codeUserVec(nextUserVec, maxUsersVec);
            
            % User departure due to finishing service 
            trm(idx, nextIdx) = userVec(ii)* ...
                (rateMat(ii,idx)/vQualMat(ii,idx))/averageVideoSizeVec(ii);
            
            % User balking due to poor video quality
            trm(idx, nextIdx) = trm(idx, nextIdx) + userVec(ii)* ...
                (vQualMat(ii,idx) < thresVec(ii)) ...
                * (gammaVec(ii) / vQualMat(ii,idx));  % may not be correct
        end
    end    
    trm(idx,idx) = -sum(trm(idx,:)); % sum of the row why?
end
