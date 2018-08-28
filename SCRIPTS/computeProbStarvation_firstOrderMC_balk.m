function [probStarVec]= computeProbStarvation_firstOrderMC_balk(trm, ...
    badStatesIdx, absorbingState, vQualVec, throughputVec, prefetch, bVec, ...
    slope, secsPerSeg, lVec, maxUsersVec, classId, probFinishing)
% Compute prob. starvation corresponding to each entry state.
% We account for first order effect, i.e., we consider the possibility of 
% transiting from a bad state to a good state before buffer depletes to zero.

% INPUT:
% trm: tagged MC trm of some class
% badStatesIdx: badStatesIdx(ii) is 1 if ii is a bad state
% absorbingState: index of absorbing state
% vQualVec: vector of avg. video bit rates
% rateVec: vector of avg. channel rates
% prefetch: prefetch threshold (number of segments)
% bVec: vector of buffer thresholds
% slope
% secsPerSeg: seconds of video in each segment
% lVec: vector of video bit-rates
% maxUsersVec:
% classId: 
% probFinishing:

% OUTPUT: 
% probStarVec: prob. starvation 


tpm = trm2Embeddedtpm(trm);
allStates = 1:length(badStatesIdx);
badStates = allStates(badStatesIdx == 1);
goodStates = allStates(badStatesIdx == 0);
probStarVec = zeros(length(allStates),1);
probStarVec(absorbingState) = 0;

wBadNeighbours = [];
for ii = goodStates,
    if sum(tpm(ii,badStates) > 0)
        wBadNeighbours = [wBadNeighbours, ii];
    end
end

wGoodNeighbours = [];
for ii = badStates,
    if sum(tpm(ii,goodStates) > 0)
        wGoodNeighbours = [wGoodNeighbours, ii];
    end
end
probStarVec(setdiff(badStates,wGoodNeighbours)) = 1;

% Find probability of avoiding starvation due to departure before
% 'prefetching'
T_B0 = zeros(size(allStates));
pB0= zeros(size(allStates));
%slope = ((bmax - bmin)/(lmax - l0));

for ii = wGoodNeighbours,
    b = prefetch*secsPerSeg;
    userVec = getUserVec(ii, maxUsersVec); % ii --> (i_1, i_2, ..., i_K)
    actualUserVec = userVec; % tagged user does not see itself in system!
    actualUserVec(classId) = userVec(classId) + 1;
    iiPlusOne = codeUserVec(actualUserVec, maxUsersVec);
    T_B0(ii) = b/((throughputVec(iiPlusOne))/min(lVec)) + ...
        b/(1 - (throughputVec(iiPlusOne))/min(lVec));
    pB0(ii) = exp(trm(ii,ii)*T_B0(ii));
end;

% Find probability of avoiding starvation due to departure before buffer
% gets empty.
pB1 = cell(size(allStates));
idVec = 1:length(lVec);
idx = length(lVec);
for ii = wBadNeighbours,
    userVec = getUserVec(ii, maxUsersVec); % ii --> (i_1, i_2, ..., i_K)
    actualUserVec = userVec; % tagged user does not see itself in system!
    actualUserVec(classId) = userVec(classId) + 1;
    iiPlusOne = codeUserVec(actualUserVec, maxUsersVec);
    if (length(bVec) > 1)
        if (vQualVec(iiPlusOne) < max(lVec))
            idx = max(idVec(vQualVec(iiPlusOne) >= lVec));
            b = bVec(idx+1)*secsPerSeg;
        else
            b = bVec(idx)*secsPerSeg;
        end
    else
        b = (slope * (vQualVec(iiPlusOne) - min(lVec)) + bVec)*secsPerSeg;
    end
    for jj = badStates,
        if(tpm(ii,jj) > 0)
            userVec = getUserVec(jj, maxUsersVec); % ii --> (i_1, i_2, ..., i_K)
            actualUserVec = userVec; % tagged user does not see itself in system!
            actualUserVec(classId) = userVec(classId) + 1;
            jjPlusOne = codeUserVec(actualUserVec, maxUsersVec);
            T_B1 = b/(1 - (throughputVec(jjPlusOne))/min(lVec));
            pB1{ii}(jj) = exp(trm(jj,jj)*T_B1);
        end;
    end;
end;

A = zeros(length(allStates));
b = zeros(length(allStates), 1);
for ii = goodStates,
    A(ii,goodStates) = -tpm(ii,goodStates);
    A(ii,ii) = 1;
end
for ii = wBadNeighbours,
    for jj = badStates,
        if(tpm(ii,jj) > 0)
            b(ii) = b(ii) + tpm(ii,jj)*pB1{ii}(jj) .* probFinishing(jj);
            b(ii) = b(ii) + tpm(ii,jj)*(1 - pB1{ii}(jj))*sum(tpm(jj,badStates).* probFinishing(badStates));
            for kk = wBadNeighbours,
                if(tpm(jj,kk) > 0)
                    C_temp = tpm(ii,jj)*(1 - pB1{ii}(jj))*tpm(jj,kk);
                    A(ii,kk) = A(ii,kk) - C_temp;
                end;
            end;
        end;
    end;
end

for ii = wGoodNeighbours,
    b(ii) = pB0(ii).* probFinishing(ii)  + ...
        (1 - pB0(ii))*sum(tpm(ii,badStates) .* probFinishing(badStates));
    A(ii,ii) = 1;
    for kk = goodStates,
        if(tpm(ii,kk) > 0)
            A(ii,kk) = A(ii,kk) - (1 - pB0(ii))*tpm(ii,kk);
        end;
    end;
end
Aye = A([goodStates, wGoodNeighbours], [goodStates, wGoodNeighbours]);
bee = b([goodStates, wGoodNeighbours]');
probStarVec([goodStates, wGoodNeighbours]') = Aye\bee;
