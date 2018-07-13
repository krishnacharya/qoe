function [avgQualVec, avgQualSwitchVec ] = computeAvgQuality_firstOrderMC_balk(trm, badStatesIdx, vQVec, rateMat, secsPerSeg, videoRateVec, probFinishing, maxUsersVec, classId)
% Code to compute analytically average video bitrate and average number of switches per class.
% INPUT:
% OUTPUT: 
tpm = trm2tpm(trm);
allStates = 1:length(badStatesIdx);
badStates = allStates(badStatesIdx == 1);
goodStates = allStates(badStatesIdx == 0);
nonAbsorbingStates = [goodStates, badStates];
numNonAbsStates = length(nonAbsorbingStates);
I = eye(numNonAbsStates);
P = tpm(nonAbsorbingStates, nonAbsorbingStates);
A_inv = zeros(length(badStatesIdx));
A_inv(nonAbsorbingStates, nonAbsorbingStates) = (I - P)^-1;
avgQualVec = zeros(1, length(allStates) - 1);
effvQVec = zeros(1,numNonAbsStates);
for ii = nonAbsorbingStates,
    count = 1;
    for jj = nonAbsorbingStates,
        userVec = getUserVec(jj, maxUsersVec); % ii --> (i_1, i_2, ..., i_K)
        actualUserVec = userVec; % tagged user does not see itself in system!
        actualUserVec(classId) = userVec(classId) + 1;
        jjPlusOne = codeUserVec(actualUserVec, maxUsersVec);
        effvQVec(count) = vQVec(jjPlusOne);
        count = count + 1;
    end;
    avgQualVec(ii) = sum( (A_inv(ii,nonAbsorbingStates) .*  ...
        probFinishing(nonAbsorbingStates) ...
        /...
        sum(A_inv(ii,nonAbsorbingStates) .* probFinishing(nonAbsorbingStates)))...
        .*effvQVec);
end;

lmin = min(videoRateVec);
lmax = max(videoRateVec);
S = zeros(length(badStatesIdx), 1);
% slope = ((bmax - bmin)/(lmax - lmin));
for ii = nonAbsorbingStates,
    userVec = getUserVec(ii, maxUsersVec); % ii --> (i_1, i_2, ..., i_K)
    actualUserVec = userVec; % tagged user does not see itself in system!
    actualUserVec(classId) = userVec(classId) + 1;
    iiPlusOne = codeUserVec(actualUserVec, maxUsersVec);
    r = rateMat(iiPlusOne);
    if ((r > lmin) && (r < lmax))
        l1 = max(videoRateVec(videoRateVec <= r));
        l2 = min(videoRateVec(videoRateVec >= r));
        if ( (l2 > l1) )
            if (r < (l1+l2)/2)
                T = secsPerSeg*(l2/r + (l2 - r)*l1/(r*(r - l1)));
            else
                T = secsPerSeg*(l1/r + (r - l1)*l2/((l2 - r)*r));
            end
            S(ii) = exp(-abs(trm(ii,ii))*T/2)/(1.0 - exp(-abs(trm(ii,ii))*T/2));
        end;
	S(ii) = S(ii)*probFinishing(ii);
    end;
end
epsilon = 1e-4;
Z = zeros(length(allStates));
for ii = nonAbsorbingStates,
    for jj = nonAbsorbingStates,
        if(tpm(ii,jj) > epsilon*epsilon)
            userVec = getUserVec(ii, maxUsersVec); % ii --> (i_1, i_2, ..., i_K)
            actualUserVec = userVec; % tagged user does not see itself in system!
            actualUserVec(classId) = userVec(classId) + 1;
            iiPlusOne = codeUserVec(actualUserVec, maxUsersVec);
            
            userVec = getUserVec(jj, maxUsersVec); % ii --> (i_1, i_2, ..., i_K)
            actualUserVec = userVec; % tagged user does not see itself in system!
            actualUserVec(classId) = userVec(classId) + 1;
            jjPlusOne = codeUserVec(actualUserVec, maxUsersVec);
            
            if(abs(vQVec(jjPlusOne) - vQVec(iiPlusOne)) > epsilon)
                Z(ii,jj) = max( sum(vQVec(jjPlusOne) <= videoRateVec ...
                    & videoRateVec <= vQVec(iiPlusOne)), ...
                    sum(vQVec(iiPlusOne) <= videoRateVec & ...
                    videoRateVec <= vQVec(jjPlusOne)));
                idx = sum(abs(vQVec(iiPlusOne) - videoRateVec) < epsilon);
                idx = idx*sum(abs(vQVec(jjPlusOne) - videoRateVec) < epsilon);
                Z(ii,jj) = Z(ii,jj) - idx;
            end
        end
    end;
    Z(ii,jj) = Z(ii,jj)*probFinishing(jj);
end

avgQualSwitchVec = zeros(1, length(allStates) - 1);
emTpm = trm2Embeddedtpm(trm);
P2 = emTpm(nonAbsorbingStates, nonAbsorbingStates);
b = S(nonAbsorbingStates,1) + (P2.*Z(nonAbsorbingStates,nonAbsorbingStates))*ones(length(nonAbsorbingStates), 1);
avgQualSwitchVec(nonAbsorbingStates) = (I - P2)^-1 * b ...
    ./  probFinishing(nonAbsorbingStates)';
