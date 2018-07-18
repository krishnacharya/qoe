function [pye, probBlocking, probFinishing, probVB, probStarvation, avgQualSwitches, ...
    avgQuality, prefetchDelay, prefetchDelayclassState] = ...
    userMC_firstOrderMC_PF_balk(arrivalRateVec, avgVideoSizeVec, ...
    channelRate, gammaVec, thresVec, weightVec, maxUsersVec, ...
    videoRateMat, prefetchVec, secsPerSegVec, bminVec, bmaxVec, unifVec)
% AUTHOR: SUDHEER
% LAST MODIFIED: 27 JUNE 2018

% This code returns analytical results for PF sharing.
% All input parameters which are vector valued are suffixed Vec and matrices
% are suffixed Mat.
% Read individual analysis file comments for details.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT PARAMETER LIST      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% arrivalRateVec: user arrival rate vectors (# of users per sec for each class)  

% avgVideoSizeVec: avg. video size (sec) for each class

% channelRate: nominal channel rate  

% gammaVec: balking rate vector (higher values corr. to  higher impatience)

% thresVec: threshold at each user gets impatient  

% weightVec: % Weight vector for max-channel rate
% A wtVec = [1, 2] and cRate of 5e6 will set max. aggregate channel rate for 
% class 1 = 5e6 and for class 2 = 10e6.

% maxUsersVec: Max. number of users per class [N_1, N_2, ...]
% N_k = max. number of users of class k that can be there in the system.

% videoRateMat: Matrix of Available bit-rates
% Row-i of videoRateMatrix gives the available bit rates for class-i.
% The vector of rates must be in ascending order.
% Append vector with zeros if the lengths are not equal

% prefetchVec: DASH parameter number of video segments to be prefetched before
% playput

% secsPerSegVec: DASH parameter, number of seconds per segment

% bminVec: DASH parameter, min. buffer threshold 

% bmaxVec: DASH parameter, max. buffer threshold 

% unifVec: Parameter which decides buffer thresholds, 0: linear, 1: thresholds uniformly 
% spaced, 2: minimum required


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  OUTPUT VARIABLES LIST      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All variables are vectors giving values per class with the exception of pye.
% pye is a vector giving stationary distribution of users.

% pye: stationary distribution of users. We list all states as a vector, where 
% each state is encoded as a scalar by function codeUserVec().
 
% probBlocking: prob. of user blocking 

% probFinishing: prop. of users who after entering service finish service

% probVB: prob. of visiting bad states (zero order prob. starvation)  

% probStarvation: prob. of starvation  (first order prob. starvation)

% avgQualSwitches: average number of quality switches 

% avgQuality: average video-bitrate 

% prefetchDelay: prefetch/startup delay 


% MAIN code for ANALYSIS
numClasses = length(maxUsersVec);
numStates = prod(maxUsersVec + 1);
epsilon = 1e-5;

% Compute steady state distribution for User Markov Chain
lminVec = zeros(1,numClasses);
lmaxVec = zeros(1,numClasses);

sim_time = clock;
filename = sprintf('analysis_PF_results_%02d-%02d-%02d-%02d-%02d.txt', sim_time(1:5));
fid = fopen(filename,'a+');

for ii=1:numClasses
    lminVec(ii) = min( removeZeros(videoRateMat(ii,:), epsilon) );
    lmaxVec(ii) = max( removeZeros(videoRateMat(ii,:), epsilon) );
end;

vQualMat = zeros(numClasses,numStates);
rateMat = zeros(numClasses,numStates);
for jj = 2:numStates,
    userVec = getUserVec(jj, maxUsersVec); % idx --> (i_1, i_2, ..., i_K)
    nr_lhs = weightVec(userVec > 0)*channelRate;
    dr_lhs = sum(userVec); % price for PF
    rateMat((userVec > 0), jj) = nr_lhs/dr_lhs;
    for ii=1:numClasses,
        if(userVec(ii) > 0)
            vQualMat(ii, jj) = max(lminVec(ii), min(lmaxVec(ii), ...
                rateMat(ii, jj)));
        end
    end;
end

trm = encodeTRM_balk(arrivalRateVec, avgVideoSizeVec, gammaVec, thresVec, ...
    vQualMat, rateMat, maxUsersVec);
pye = getStationaryDist(trm2tpm(trm));
probBlocking = computeProbBlocking(pye, maxUsersVec);

% Compute TRMs for tagged user MC
trm2 = cell(1,numClasses);
for kk=1:numClasses,
    % TRM for tagged user MC of class kk, each class has its own transition rate matrix
    trm2{kk} = encodeTaggedTRM_balk(arrivalRateVec, avgVideoSizeVec, ...
        gammaVec, thresVec, vQualMat, rateMat, maxUsersVec, kk);
end;

% Code the badStates as +1, if the state is not possible, encode it as -1,
% goodStates are 0, absorbing states are 2 and 3
badStates = encodeBadStates_balk(vQualMat, rateMat, maxUsersVec);
absorbingState = codeUserVec(maxUsersVec, maxUsersVec);% no more of any class users possible?
badStates(:,absorbingState) = 2; % Finish service
badStates(:,absorbingState+1) = 3; % Abort service

% Compute probability of finishing service for each class, all initial
% states
probFinishing = zeros(1,numClasses);
probFinishingMat = zeros(numClasses, numStates);
pyeTimesprobFinishingMat = zeros(numClasses, numStates);
for kk=1:numClasses,
    % TRM for tagged user MC of class kk
    listStates = 1:length(badStates(kk,:)); % this is just the number of states (N1+1) * (N2+1) etc
    validStates = listStates(badStates(kk,:) == 0 | badStates(kk,:) == 1); % valid state can be good or bad.
    probFinishingMat(kk, validStates) = computeProbFinishing(trm2{kk}, badStates(kk,:));
    temp = probFinishingMat(kk, :);
    pyeTimesprobFinishingMat(kk, :) = pye' .* temp; %pye is 21x1 in shape 
    probFinishing(kk) = sum(pyeTimesprobFinishingMat(kk, :));
    probBlocking(kk) = probBlocking(kk) / (probFinishing(kk) + probBlocking(kk));
end;


% Compute starvation probability [UPPER BOUND]
probVB = zeros(1,numClasses);
for kk=1:numClasses,
    tpm = trm2Embeddedtpm(trm2{kk});
    probVBVec = zeros(1,numStates);
    probVBVec(badStates(kk,:) == 1) = 1;
    [Ps, goodStates] = ...
        computeProbStarvation_zeroOrderMC_balk(tpm, badStates(kk,:), ...
        probFinishingMat(kk, :));
    probVBVec(goodStates) = Ps' ./ probFinishingMat(kk,goodStates) ;
    probVB(kk) = sum(pyeTimesprobFinishingMat(kk,:).*probVBVec) / ...
        sum(pyeTimesprobFinishingMat(kk,:));
end;

% Compute starvation probability [MORE EXACT BOUND]
bVec = cell(1,numClasses);
probStarvation = zeros(1,numClasses);
for kk=1:numClasses,
    bVec{kk} = computeBufferSpacing(videoRateMat(kk,:), bminVec(kk), bmaxVec(kk), unifVec(kk));
    if(unifVec(kk) > 2)
        bVec{kk} = bminVec(kk);
    end
    slope = (bmaxVec(kk) * bminVec)/(lmaxVec(kk) - lminVec(kk));
    probStarVec = computeProbStarvation_firstOrderMC_balk(trm2{kk}, ...
        badStates(kk,:), absorbingState, vQualMat(kk,:), rateMat(kk,:), ...
        prefetchVec(kk), bVec{kk}, slope, secsPerSegVec(kk), ...
        videoRateMat(kk,:), maxUsersVec, kk, probFinishingMat(kk, :));
    probStarVec(badStates(kk,:) <= 1)  = probStarVec(badStates(kk,:) <= 1) ./ ...
        probFinishingMat(kk,badStates(kk,:) <= 1)';
    probStarvation(kk) = sum(pyeTimesprobFinishingMat(kk,badStates(kk,:) <= 1)...
        .*probStarVec(badStates(kk,:) <= 1)') / sum(pyeTimesprobFinishingMat(kk,:));
end;

% Compute average prefetchDelay
prefetchDelay = zeros(1,numClasses);
prefetchDelayclassState = zeros(numClasses, numStates);
for kk=1:numClasses,
    effRateVec = zeros(1, numStates);
    for ii = 1:numStates,
        userVec = getUserVec(ii, maxUsersVec); % ii --> (i_1, i_2, ..., i_K)
        if (userVec(kk) < maxUsersVec(kk))
            actualUserVec = userVec; % tagged user does not see itself in system!
            actualUserVec(kk) = userVec(kk) + 1;
            iiPlusOne = codeUserVec(actualUserVec, maxUsersVec);
            effRateVec(ii) = rateMat(kk,iiPlusOne);% r(i + ej)
            prefetchDelayclassState(kk, ii) = prefetchVec(kk) * lminVec(kk) / rateMat(kk, iiPlusOne);
        end;
    end;
    
    prefetchDelay(kk) = sum( (prefetchVec(kk)*secsPerSegVec(kk)* ...
        lminVec(kk) ./ effRateVec(effRateVec > epsilon)) ...
        ...
        .*  pyeTimesprobFinishingMat(kk,(effRateVec > epsilon))...
        / sum(pyeTimesprobFinishingMat(kk, effRateVec > epsilon)) );
end;


% Compute average video quality and quality switches
avgQualSwitches = zeros(1,numClasses);
avgQuality = zeros(1,numClasses);
for kk=1:numClasses,
    [avgQualVec, avgQualSwitchVec ]= computeAvgQuality_firstOrderMC_balk(...
        trm2{kk}, badStates(kk,:), vQualMat(kk,:), rateMat(kk,:), ...
        secsPerSegVec(kk), videoRateMat(kk,:), probFinishingMat(kk,:), ...
        maxUsersVec, kk);    
    avgQualSwitches(kk) = sum(pyeTimesprobFinishingMat(kk,:) .* ...
        avgQualSwitchVec) / (sum(pyeTimesprobFinishingMat(kk,:)));
    avgQuality(kk) = sum(pyeTimesprobFinishingMat(kk,:) .* ...
        avgQualVec)/(sum(pyeTimesprobFinishingMat(kk,:)));
end



fprintf(fid,'\nAnalysis steady state distr: ');
printVec(fid, pye', length(pye));

for jj = 1:length(maxUsersVec),
    fprintf(fid,'\nAnalysis FINAL Res: class %d ', jj);
    fprintf(fid,'S: %.5f, B: %.5f, Q: %.3f, NQS: %.3f, RQS: %.3f, D: %.3f, PD: %.3f, PF: %.4f', ...       
        probStarvation(jj), probVB(jj), avgQuality(jj), avgQualSwitches(jj), ...
        avgQualSwitches(jj)/avgVideoSizeVec(jj), prefetchDelay(jj), ...
        probBlocking(jj), (sum(pyeTimesprobFinishingMat(jj,:))) / ...
        (1.0 - probBlocking(jj)) );
end;
fprintf(fid,'\n\n');
fclose(fid);
