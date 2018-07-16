function [pi, pi_2, user,   probStarvClass, probVBClass, probFinishClass, probDropClass, avgQualityClass, avgQualitySwitchesClass, ...
    avgPrefetchTimeClass, avgDownloadTimeClass, avgVideoDurationClass, numUsers, simTime] = ...
    script_Imen_SP_MC_PF_balk(arrivalRateVec, avgVideoSizeVec, ...
    channelRate, gammaVec, minRateThresVec, weightVec, maxUsersVec, videorateMatrix, ...
    unifVec, avgUsersSim)
% balking added to Proportional fair sharing
% unif for buffer spacing
% 0: linear, 1: unif. spacing, 2: minimum reqd., else: ctx. video quality

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT parameters description
% arrivalRateVec: User arrival rate vector (per class in arrivals/sec)
% avgVideoSizeVec: Average video sizei vector (per class in seconds)
% channelRate: Nominal channel rate
% gammaVec: Vector of multipliers for balking, higher values give higher balking rate
% minRateThresVec: Vector or minimum rate required for user satisfaction, users MAY leave if rate received is less than this
% weightVec = 1: % Weight vector for max-channel rate
	% A wtVec = [1, 2] and cRate of 5e6 will set max. aggregate channel rate for class 1 = 5e6 and for class 2 = 10e6.
% maxUsersVec: Max. number of users per class
% videoRateMatrix: Vector of Available bit-rates
	% Row-i of videoRateMatrix gives the available bit rates for class-i.
	% The vector of rates must be in ascending order.
	% Append vector with zeros if the lengths are not equal
%i unifVec: Parameter which decides buffer thresholds, 0: linear, 1: thresholds uniformly spaced, 2: minimum required
% avgUsersSim: The simulation will simulate (on average) avgUsersSim users entering the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT variables description
% pi: Distribution of users in system as seen by incoming users 
% pi_2: Distribution of users in system per simulation slot
% probStarvClass: prob. of starvation per class
% probVBClass: prob. of user visiting a `bad' state, i.e., state with available rate less than l_{min}
% probFinishClass: prob. of user finishing service
% probDropClass: prob. that user was dropped (considering only users who finished service)
% avgQualityClass: avg. video bit-rate per class 
% avgQualitySwitchesClass: avg. number of quality switches per class
% avgPrefetchTimeClass: avg. prefetch delay per class
% avgDownloadTimeClass: avg. time to finish video download
% avgVideoDurationClass: avg. duration of video
% numUsers: number of users considered for computing simulation results
% simTime: simulation duration (duration of simulation is NOT the actual running time of the program)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% rng('shuffle');
rand('state', 42) ;%just seeding with 42
paramsDASH(1:length(maxUsersVec)) = struct;
minRateVec = zeros(1,length(maxUsersVec));
%% DASH parameter settings
for jj=1:length(maxUsersVec)
    paramsDASH(jj).lambda = 30 ; % frame rate in fps
    paramsDASH(jj).qs= 1; % buffer cache/ prefetching threshold for startup in (segments)
    paramsDASH(jj).segment = 60; % segment size in (frames)
    paramsDASH(jj).l= removeZeros(videorateMatrix(jj,:), 1e-1); %[0.2 0.3 0.48 0.75 1.2 1.85 2.85 4.3 5.3]*1e6; % list of quality levels
    paramsDASH(jj).Bmin=4;
    paramsDASH(jj).Bc=10;
    paramsDASH(jj).segPerSec = paramsDASH(jj).lambda/paramsDASH(jj).segment;
    minRateVec(jj) = min(paramsDASH(jj).l);
end

badStates = zeros(length(maxUsersVec), prod(maxUsersVec + 1));
for jj = 1:prod(maxUsersVec + 1),
    userVec = getUserVec(jj, maxUsersVec); % idx --> (i_1, i_2, ..., i_K)
    nr_lhs = weightVec(userVec > 0)*channelRate;
    dr_lhs = sum(userVec);
    lhs = nr_lhs/dr_lhs; % rate for PF
    badStates((userVec > 0),jj) = lhs < minRateVec(userVec > 0);
end

%% Simulation Parameters
%simSlotsPerSec = max(10, round(10*min(arrivalRateVec)));%change this to max or sum of all arr
simSlotsPerSec = ceil(10 * sum(arrivalRateVec));
simTime = avgUsersSim / sum(arrivalRateVec); % simulation duration in seconds
simulationDuration = ceil(simTime*simSlotsPerSec); % in sim. slots
startTime = floor(0.1*simulationDuration); % discard initial 10 % samples
userThreshold = 2*sum(maxUsersVec); % discard userThreshold users
segPerSlotVec = zeros(1,length(maxUsersVec));

bThres = cell(1,length(maxUsersVec));
idThres = cell(1,length(maxUsersVec));
for jj=1:length(maxUsersVec)
    segPerSlotVec(jj) = paramsDASH(jj).segPerSec/simSlotsPerSec;
    %[IMP: Need to set video quality at the same scale as channel rate]
    paramsDASH(jj).l=paramsDASH(jj).l/simSlotsPerSec; % quality level in bp sim. slot
    paramsDASH(jj).rate = (paramsDASH(jj).l(end)-paramsDASH(jj).l(1))/...
        (paramsDASH(jj).Bc-paramsDASH(jj).Bmin);
    paramsDASH(jj).c = paramsDASH(jj).l(end) - paramsDASH(jj).rate*paramsDASH(jj).Bc;
    bThres{jj} = computeBufferSpacing(paramsDASH(jj).l, paramsDASH(jj).Bmin, paramsDASH(jj).Bc, unifVec(jj));
    idThres{jj} = 1:length(bThres{jj});
end


%% the arrival events follow the poisson distribution
% interarrival time in sim. slots
nextArrival = 1 + round(simSlotsPerSec*exprnd(1/sum(arrivalRateVec)));%check this no need for reciprocal

nextBalk = 0;
effBalkRateVec = zeros(size(gammaVec));
%% structure of user
% states :
% 1) "requesting" : the user must decide on the quality
% 2) "streaming" : the user is juste downloading the curentSegment
% 3) "outOfSystem" : the user has finished the video and is no more in the
% loop of active users

% The first user(s) add one user of each class
userVec = zeros(1, length(maxUsersVec));
idx = 1;
for idx= 1:length(maxUsersVec)
    u_= User;
    u_.id = idx;
    u_.class = idx;
    u_.videoLength=exprnd(avgVideoSizeVec(idx)); % in seconds
    u_.totalNbrOfSegments=max(1, round(u_.videoLength*paramsDASH(idx).segPerSec));
    u_.usersInSystem = userVec;
    userVec(idx) = userVec(idx) + 1;
    user(idx) = u_;% So in our user arrray we now have 1...k indices and each has a user of that class index 
end;

%C_0: Total capacity available per sim. slot
%C: Capacity available per sim. slot per user
C_0 = (channelRate/simSlotsPerSec); % bits available per sim slot
minRateThresVec = minRateThresVec/simSlotsPerSec;%if users of class i receive less than min..(i) bits then they may leave
dropCnt = 0;

s_idx = 1; % id of first user still in system
numDroppedVec = zeros(1, length(maxUsersVec));
countDroppedVec = zeros(1, length(maxUsersVec));
numUsersMatrix = zeros(simulationDuration, length(maxUsersVec));
fprintf('\n');
for s = 1:simulationDuration,  % in simulation slots
    %     if(mod(s,0.1*simulationDuration) == 0)
    %         fprintf('!');
    %     end;
    if (nextArrival == s) % Arrival occurs
        nextArrival = nextArrival + ...
            max(1, round(simSlotsPerSec * exprnd(1/sum(arrivalRateVec)))) %this again should be changed no need for reciprocal
        userType = selectUserType(arrivalRateVec);
        
        if (userVec(userType) < maxUsersVec(userType))
            idx = idx + 1;
            u_ = User;
            u_.class = userType;
            u_.id= idx;
            u_.videoLength = exprnd(avgVideoSizeVec(userType));
            u_.totalNbrOfSegments = max(1, round(u_.videoLength * paramsDASH(userType).segPerSec)); % Video length in segments
            u_.usersInSystem = userVec;
            u_.entryTime = s;
            userVec(userType) = userVec(userType) + 1;
            
            % User balking rate changed:
            % Change balk time as num. users changed
            dr_lhs = sum(userVec);
            rateVec = zeros(1,length(maxUsersVec));
            rateVec(userVec > 0) = weightVec(userVec > 0)*C_0;
            rateVec = rateVec/dr_lhs;
            
            effBalkRateVec = gammaVec .* userVec ./ (rateVec .* simSlotsPerSec) ...
                            .* (rateVec < minRateThresVec);
            effBalkRate = sum(effBalkRateVec);
            if(effBalkRate > 0)
                nextBalk = s + max(1, round(simSlotsPerSec*exprnd(1/effBalkRate)));
            else
                nextBalk = simulationDuration + 1;
            end;
        else
            numDroppedVec(userType) = numDroppedVec(userType) + 1;
            if(length(user) + sum(numDroppedVec) >= userThreshold)
                u_ = User;
                u_.class = userType;
                u_.usersInSystem = userVec;
                u_.entryTime = s;
                dropCnt = dropCnt + 1;
                droppedUsers(dropCnt) = u_;
                countDroppedVec(userType) = countDroppedVec(userType) + 1;
            end;
        end
    end
    
    if (nextBalk == s)
        userType = selectUserType(effBalkRateVec);
        if(userVec(userType) > 0 && userType > 0)
            cntResSampling = 0;
            for i = s_idx:length(user)
                if ((user(i).class == userType) && (user(i).state < 2))
                    if(rand < 1/(userVec(userType) - cntResSampling))
                        user(i).state = 3;
                        userVec(userType) = userVec(userType) - 1;
                        break;
                    end
                    cntResSampling = cntResSampling+1;
                end;
            end;
            
            % Set new balk time
            dr_lhs = sum(userVec);
            rateVec = zeros(1, length(maxUsersVec));
            rateVec(userVec > 0) = weightVec(userVec > 0) * C_0;% this will have the rates for all the classes at the current state
            rateVec = rateVec / dr_lhs;
            
            effBalkRateVec = gammaVec .* userVec ./ (rateVec .* simSlotsPerSec) ...
                            .* (rateVec < minRateThresVec);
            effBalkRate = sum(effBalkRateVec);
            if(effBalkRate > 0)
                nextBalk = s + max(1, round(simSlotsPerSec*exprnd(1/effBalkRate)));
            else
                nextBalk = simulationDuration + 1;
            end;
        end;
    end;
    
    
    if (sum(userVec) > 0) % Whenever there is a user in the system
        s_flag = 0; % flag which indicates s_idx changed.
        dr_lhs = sum(userVec);
        rateVec = zeros(1,length(maxUsersVec));
        rateVec(userVec > 0) = weightVec(userVec > 0) * C_0;
        rateVec = rateVec/dr_lhs;
        depFlag = 0;
        prevUserVec = userVec;
        for i = 1:length(user)
            if  (user(i).state < 2)
                if (0 == s_flag)%worthless check
                    s_idx = i;
                    s_flag = 1;
                end;
                if ( badStates(user(i).class, codeUserVec(userVec, maxUsersVec) ) == 1)
                    user(i).visitsBadState = 1;
                end
                [user(i), userVec(user(i).class)] = DASH_KR(paramsDASH(user(i).class), segPerSlotVec(user(i).class), user(i), s, ...
                    rateVec(user(i).class), userVec(user(i).class), bThres{user(i).class}, idThres{user(i).class});                
                if(prevUserVec(user(i).class) < userVec(user(i).class)) %check this
                    depFlag = 1;
                end;
            end
        end
    end
    numUsersMatrix(s,:) = userVec;
    
    % User balking rate changed:
    if ( (s > 1) && (depFlag == 1) )        
        % Change balk time as num. users changed
        dr_lhs = sum(userVec);
        rateVec = zeros(1,length(maxUsersVec));
        rateVec(userVec > 0) = weightVec(userVec > 0)*C_0;
        rateVec = rateVec/dr_lhs;
        
        effBalkRateVec = gammaVec .* userVec ./ (rateVec .* simSlotsPerSec) ...
                            .* (rateVec < minRateThresVec);
        effBalkRate = sum(effBalkRateVec);
        if(effBalkRate > 0)
            nextBalk = s + max(1, round(simSlotsPerSec*exprnd(1/effBalkRate)));
        else
            nextBalk = simulationDuration + 1;
        end;
    end;
    
end


%% Compute empirical distribution of number of users seen by incoming users
pi = zeros(1,prod(maxUsersVec + 1));
% We ignore the first userThreshold users
startUser = userThreshold - sum(numDroppedVec) + sum(countDroppedVec);
for ii=startUser:length(user)
    userVec = user(ii).usersInSystem;
    pi(codeUserVec(userVec, maxUsersVec)) = pi(codeUserVec(userVec, maxUsersVec)) + 1;
end

if(dropCnt > 0)
    for ii=1:length(droppedUsers)
        userVec = droppedUsers(ii).usersInSystem;
        pi(codeUserVec(userVec, maxUsersVec)) = pi(codeUserVec(userVec, maxUsersVec)) + 1;
    end
end

if(sum(pi) > 0)
    pi = pi/(sum(pi));
end

%% Compute empirical distribution of number of users (per slot granularity)
pi_2 = zeros(1,prod(maxUsersVec + 1));
for s=startTime:simulationDuration,
    userVec = numUsersMatrix(s,:);
    pi_2(codeUserVec(userVec, maxUsersVec)) = pi_2(codeUserVec(userVec, maxUsersVec)) + 1;
end
if(sum(pi_2) > 0)
    pi_2 = pi_2/(sum(pi_2));
end

%% Compute probability of starvation, # quality switches, avg quality, prefetching time
countFinishVec = zeros(1,length(maxUsersVec));
countBalkVec = zeros(1,length(maxUsersVec));
starvCntVec = zeros(1,length(maxUsersVec));
vBCntVec = zeros(1,length(maxUsersVec));
avgQualitySwitchesClass = zeros(1,length(maxUsersVec));
avgQualityClass = zeros(1,length(maxUsersVec));
avgPrefetchTimeClass = zeros(1,length(maxUsersVec));
avgDownloadTimeClass = zeros(1,length(maxUsersVec));
avgVideoDurationClass = zeros(1,length(maxUsersVec));


for i=startUser:length(user),
    if(user(i).state == 3)
        idx = user(i).class;
        countBalkVec(idx) = countBalkVec(idx) + 1;
    end;
    if(user(i).state == 2)
        idx = user(i).class;
        countFinishVec(idx) = countFinishVec(idx) + 1;
        
        avgQualitySwitchesClass(idx) = avgQualitySwitchesClass(idx) + user(i).qualitySwitches - 1;
        
        avgPrefetchTimeClass(idx) = avgPrefetchTimeClass(idx) + ...
            (user(i).prefetchTime - user(i).entryTime)/simSlotsPerSec;
        
        avgDownloadTimeClass(idx) = avgDownloadTimeClass(idx) + ...
            (user(i).exitTime - user(i).entryTime)/simSlotsPerSec;
        
        avgVideoDurationClass(idx) = avgVideoDurationClass(idx) + ...
            user(i).totalNbrOfSegments/paramsDASH(idx).segPerSec;
        
        avgQualityClass(idx) = avgQualityClass(idx) + user(i).avgQuality;
        
        if(user(i).starv > 0)
            starvCntVec(idx) = starvCntVec(idx) + 1;
        end
        
        if(user(i).visitsBadState > 0)
            vBCntVec(idx) = vBCntVec(idx) + 1;
        end
        
    end
end

if(sum(countFinishVec == 0) == 0)
    avgQualitySwitchesClass = avgQualitySwitchesClass./countFinishVec;
    avgPrefetchTimeClass = avgPrefetchTimeClass./countFinishVec;
    avgDownloadTimeClass = avgDownloadTimeClass./countFinishVec;
    avgVideoDurationClass = avgVideoDurationClass./countFinishVec;
    avgQualityClass = avgQualityClass./countFinishVec; %in bit per sim. slot
    avgQualityClass = avgQualityClass*simSlotsPerSec; % in bps
    probStarvClass = starvCntVec./countFinishVec;
    probVBClass = vBCntVec./countFinishVec;
    probFinishClass = countFinishVec ./ (countFinishVec + countBalkVec);
    probDropClass = countDroppedVec ./ (countFinishVec + countDroppedVec);
else
    minusOneVec = -1*ones(1,length(maxUsersVec));
    avgQualitySwitchesClass = minusOneVec;
    avgPrefetchTimeClass = minusOneVec;
    avgDownloadTimeClass = minusOneVec;
    avgVideoDurationClass = minusOneVec;
    avgQualityClass = minusOneVec;
    probStarvClass = minusOneVec;
    probVBClass = minusOneVec;
    probFinishClass = minusOneVec;
    probDropClass = minusOneVec;
end
numUsers = sum(countFinishVec);
% save results.mat; %For error checking
