function [pi, pi_2, user, avgTimeinSystem, probStarvClass, probVBClass, probFinishClass, probDropClass, avgQualityClass, avgQualitySwitchesClass, avgPrefetchTimeClass, AvgPrefetchTimeij, FreqMatrix, avgDownloadTimeClass, avgVideoDurationClass, numUsers, simTime, avgSquaredQualDiff] = ...
    simscript2(arrivalRateVec, prefetchVec, avgVideoSizeVec, secsPerSegVec, gammaVec, minRateThresVec, ...
    maxUsersVec, videorateMatrix, unifVec, bminVec, bmaxVec, avgUsersSim, throughputVec)
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
%picking relevant variable from the Map
% arrivalRateVec = M('lambda');
% prefetchVec = M('prefetchVec');
% avgVideoSizeVec = M('avSizeVec');
% secsPerSegVec = M('secsPerSegVec');
% %channelRate 
% gammaVec = M('gammaVec'); 
% minRateThresVec = M('minRateThresVec');
% %weightVec
% maxUsersVec = M('maxUsersVec');
% videorateMatrix = M('videoRateMatrix');
% unifVec = M('unifVec'); 
% bminVec = M('bminVec');
% bmaxVec = M('bmaxVec');
% avgUsersSim = M('avgUsersSim');
% throughputVec = M('throughputVec');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state', 42) ;%just seeding with 42
numberOfClasses = length(maxUsersVec);
numberOfStates = prod(maxUsersVec + 1);
paramsDASH(1 : numberOfClasses) = struct;
minRateVec = zeros(1, numberOfClasses);
%% DASH parameter settings
for jj = 1: numberOfClasses
    paramsDASH(jj).qs = prefetchVec(jj); % buffer cache/ prefetching threshold for startup in (segments)
    paramsDASH(jj).l = removeZeros(videorateMatrix(jj,:), 1e-1);
    paramsDASH(jj).Bmin = bminVec(jj);
    paramsDASH(jj).Bc = bmaxVec(jj);
    paramsDASH(jj).segPerSec = 1.0 / secsPerSegVec(jj)
    minRateVec(jj) = min(paramsDASH(jj).l);
end

badStates = zeros(numberOfClasses, numberOfStates);%flag of 1 if it is a bad state
for jj = 2 : numberOfStates
    userVec = getUserVec(jj, maxUsersVec); % idx --> (i_1, i_2, ..., i_K)
    badStates((userVec > 0), jj) = throughputVec(userVec) < minRateVec(userVec > 0);  
end
%% Simulation Parameters
simSlotsPerSec = 10 * ceil(30 * sum(arrivalRateVec));
secPerSlot = (1.0 / simSlotsPerSec);
simTime = avgUsersSim / sum(arrivalRateVec); % simulation duration in seconds
simulationDuration = ceil(simTime * simSlotsPerSec); % in sim. slots
startTime = floor(0.1 * simulationDuration); % discard initial 10 % samples
userThreshold = 2 * sum(maxUsersVec);
segPerSlotVec = zeros(1, numberOfClasses);

bThres = cell(1, numberOfClasses);
idThres = cell(1, numberOfClasses);
for jj = 1 : numberOfClasses
    segPerSlotVec(jj) = paramsDASH(jj).segPerSec / simSlotsPerSec;
    paramsDASH(jj).l = paramsDASH(jj).l / simSlotsPerSec; % quality level in bits per sim. slot
    paramsDASH(jj).rate = (paramsDASH(jj).l(end) - paramsDASH(jj).l(1)) / (paramsDASH(jj).Bc-paramsDASH(jj).Bmin);%slope for the graph
    paramsDASH(jj).c = paramsDASH(jj).l(end) - paramsDASH(jj).rate * paramsDASH(jj).Bc;
    bThres{jj} = computeBufferSpacing(paramsDASH(jj).l, paramsDASH(jj).Bmin, paramsDASH(jj).Bc, unifVec(jj));
    idThres{jj} = 1 : length(bThres{jj});
end
%% the arrival events follow the poisson distribution
% interarrival time in sim. slots
nextArrival = 1 + round(simSlotsPerSec * exprnd(1 / sum(arrivalRateVec)));
nextBalk = 0;
effBalkRateVec = zeros(size(gammaVec));
%% structure of user
% states :
% 0) "requesting" : the user must decide on the quality
% 1) "streaming" : the user is juste downloading the curentSegment
% 2) "outOfSystem" : the user has finished the video and is no more in the
% loop of active users

% The first user(s) add one user of each class
userVec = zeros(1, numberOfClasses);
idx = 1;
for idx = 1:numberOfClasses
    u_= User;
    u_.id = idx;
    u_.class = idx;
    u_.videoLength = exprnd(avgVideoSizeVec(idx)); % in seconds
    u_.totalNbrOfSegments = max(1, round(u_.videoLength * paramsDASH(idx).segPerSec));
    u_.usersInSystem = userVec;
    userVec(idx) = userVec(idx) + 1;
    user(idx) = u_;% So in our user arrray we now have 1...k indices and each has a user of that class index 
end;

%C_0: Total capacity available per sim. slot
%C: Capacity available per sim. slot per user
%C_0 = (channelRate / simSlotsPerSec); % bits available per sim slot
minRateThresVec = minRateThresVec / simSlotsPerSec;%if users of class i receive less than min..(i) bits then they may leave
dropCnt = 0;

s_idx = 1; % id of first user still in system
numDroppedVec = zeros(1, numberOfClasses);
countDroppedVec = zeros(1, numberOfClasses);
numUsersMatrix = zeros(simulationDuration, numberOfClasses);
AvgPrefetchTimeij = zeros(numberOfClasses, numberOfStates);% average prefetch time for class i in state j
FreqMatrix = zeros(numberOfClasses, numberOfStates);
%prolongationfactor = arrivalRateVec * 100;
fprintf('\n');
for s = 1 : simulationDuration,  % in simulation slots
    if (nextArrival == s) % Arrival occurs
        nextArrival = nextArrival + max(1, round(simSlotsPerSec * exprnd(1 / sum(arrivalRateVec))))
        userType = selectUserType(arrivalRateVec);        
        if (userVec(userType) < maxUsersVec(userType))
            idx = idx + 1;
            u_ = User;
            u_.class = userType;
            u_.id = idx;
            u_.videoLength = exprnd(avgVideoSizeVec(userType));
            u_.totalNbrOfSegments = max(1, round(u_.videoLength * paramsDASH(userType).segPerSec)); % Video length in segments
            u_.usersInSystem = userVec;
            u_.entryTime = s;
            user(idx) =  u_; %adds this user to our array
            userVec(userType) = userVec(userType) + 1;
            
            %User balking rate changed:
            %Change balk time as num. users changed
            %dr_lhs = sum(userVec);
            %rateVec = zeros(1, length(maxUsersVec));
            %rateVec(userVec > 0) = weightVec(userVec > 0) * C_0; % C_0 is in bits per slot
            %rateVec = rateVec / dr_lhs;
            %rateVec = getNewRateVec(userVec, weightVec, C_0, maxUsersVec);
            
            % effBalkRateVec = gammaVec .* userVec ./ (rateVec * simSlotsPerSec) .* (rateVec < minRateThresVec);% only those with less than 				minThres seem to balk
            % effBalkRate = sum(effBalkRateVec);
            % if(effBalkRate > 0)
            %     nextBalk = s + max(1, round(simSlotsPerSec * exprnd( 1 / effBalkRate)));
            % else
            %     nextBalk = simulationDuration + 1;
            % end;            
        else
            numDroppedVec(userType) = numDroppedVec(userType) + 1;
            if(length(user) + sum(numDroppedVec) >= userThreshold)
                u_ = User;
                u_.class = userType;
                u_.usersInSystem = userVec;
                u_.entryTime = s;
                dropCnt = dropCnt + 1;
                droppedUsers(dropCnt) = u_;% stores dropped users, note they never entered the users array
                countDroppedVec(userType) = countDroppedVec(userType) + 1;
            end
        end
    end    
    % if (nextBalk == s)
    %     userType = selectUserType(effBalkRateVec);
    %     if(userVec(userType) > 0 && userType > 0)
    %         cntResSampling = 0;
    %         for i = s_idx:length(user)
    %             if ((user(i).class == userType) && (user(i).state < 2))
    %                 if(rand < 1/(userVec(userType) - cntResSampling))
    %                     user(i).state = 3;
    %                     userVec(userType) = userVec(userType) - 1;
    %                     break;
    %                 end
    %                 cntResSampling = cntResSampling+1;
    %             end;
    %         end;
           
    %         rateVec = getNewRateVec(userVec, weightVec, C_0, maxUsersVec);            
    %         effBalkRateVec = gammaVec .* userVec ./ (rateVec .* simSlotsPerSec) ...
    %                         .* (rateVec < minRateThresVec);
    %         effBalkRate = sum(effBalkRateVec);
    %         if(effBalkRate > 0)
    %             nextBalk = s + max(1, round(simSlotsPerSec*exprnd(1/effBalkRate)));
    %         else
    %             nextBalk = simulationDuration + 1;
    %         end;
    %     end;
    % end;   
    if (sum(userVec) > 0) % Whenever there are users in the system
        s_flag = 0; % flag which indicates s_idx changed.
        %rateVec = getNewRateVec(userVec, weightVec, C_0, maxUsersVec);
        depFlag = 0;
        prevUserVec = userVec;
        ScalarSystemState = codeUserVec(userVec, maxUsersVec);
        userVecUpdate = zeros(1, numberOfClasses);
        for i = 1 : length(user)
            if(user(i).state < 2)
                if (badStates(user(i).class, codeUserVec(userVec, maxUsersVec)) == 1)
                    user(i).visitsBadState = 1;
                end
                %printf('%d classOfUser \n',user(i).class);
                %to ensure it is bits per sim slot, we divide bps/ssps
                [user(i), update] = DASH_KR(paramsDASH(user(i).class), segPerSlotVec(user(i).class), user(i), s, ...
                    throughputVec(userVec) / simSlotsPerSec, bThres{user(i).class}, idThres{user(i).class}, ScalarSystemState);
                            
                userVecUpdate(user(i).class) = userVecUpdate(user(i).class) + update;
                
                if(user(i).updateDelayMatrixFlag)
                    scalarState = codeUserVec(user(i).usersInSystem,  maxUsersVec);                    
                	%AvgPrefetchTimeij(user(i).class, user(i).stateAtPrefetchStart) = AvgPrefetchTimeij(user(i).class, user	(i).stateAtPrefetchStart) + user(i).prefetchTime / simSlotsPerSec;
                    AvgPrefetchTimeij(user(i).class, scalarState) = AvgPrefetchTimeij(user(i).class, scalarState) + user(i).prefetchTime / simSlotsPerSec;
                	%FreqMatrix(user(i).class, user(i).stateAtPrefetchStart) =  FreqMatrix(user(i).class, user(i).stateAtPrefetchStart) + 1;
                    FreqMatrix(user(i).class, scalarState) =  FreqMatrix(user(i).class, scalarState) + 1;
                	user(i).updateDelayMatrixFlag = 0;
                end
            end
        end
        userVec = userVec + userVecUpdate;%update the number of users only at the end of the slot, cause they run in parallel in each slot.
        if (sum(userVecUpdate) < 0)
            depFlag = 1;
        end
    end
    numUsersMatrix(s,:) = userVec;
    
    % User balking rate changed:
    % if ( (s > 1) && (depFlag == 1) )
    %     rateVec = getNewRateVec(userVec, weightVec, C_0, maxUsersVec);
        
    %     effBalkRateVec = gammaVec .* userVec ./ (rateVec .* simSlotsPerSec) ...
    %                         .* (rateVec < minRateThresVec);
    %     effBalkRate = sum(effBalkRateVec);
    %     if(effBalkRate > 0)
    %         nextBalk = s + max(1, round(simSlotsPerSec*exprnd(1/effBalkRate)));
    %     else
    %         nextBalk = simulationDuration + 1;
    %     end;
    % end;
end
AvgPrefetchTimeij = AvgPrefetchTimeij ./ FreqMatrix;
%% Compute empirical distribution of number of users seen by incoming users
pi = zeros(1, numberOfStates);
% We ignore the first userThreshold users
startUser = userThreshold - sum(numDroppedVec) + sum(countDroppedVec);%this has to be altered
for ii = startUser:length(user)
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
pi_2 = zeros(1, numberOfStates);
for s = startTime:simulationDuration,
    userVec = numUsersMatrix(s,:);
    pi_2(codeUserVec(userVec, maxUsersVec)) = pi_2(codeUserVec(userVec, maxUsersVec)) + 1;
end
if(sum(pi_2) > 0)
    pi_2 = pi_2/(sum(pi_2));
end

%% Compute probability of starvation, # quality switches, avg quality, prefetching time
countFinishVec = zeros(1, numberOfClasses);
countBalkVec = zeros(1, numberOfClasses);
starvCntVec = zeros(1, numberOfClasses);
vBCntVec = zeros(1, numberOfClasses);
avgQualitySwitchesClass = zeros(1, numberOfClasses);
avgQualityClass = zeros(1, numberOfClasses);
avgPrefetchTimeClass = zeros(1, numberOfClasses);
avgDownloadTimeClass = zeros(1, numberOfClasses);
avgVideoDurationClass = zeros(1, numberOfClasses);
avgSquaredQualDiff = zeros(1, numberOfClasses);
for i = startUser:length(user)
    % if(user(i).state == 3)
    %     idx = user(i).class;
    %     countBalkVec(idx) = countBalkVec(idx) + 1;
    % end;            
    if(user(i).state == 2)
        idx = user(i).class;
        countFinishVec(idx) = countFinishVec(idx) + 1; 
        avgQualitySwitchesClass(idx) = avgQualitySwitchesClass(idx) + user(i).qualitySwitches - 1;        
        avgPrefetchTimeClass(idx) = avgPrefetchTimeClass(idx) + user(i).prefetchTime / simSlotsPerSec;        
        avgDownloadTimeClass(idx) = avgDownloadTimeClass(idx) + (user(i).exitTime - user(i).entryTime) / simSlotsPerSec;        
        avgVideoDurationClass(idx) = avgVideoDurationClass(idx) + user(i).videoLength;% in seconds        
        avgQualityClass(idx) = avgQualityClass(idx) + user(i).avgQuality;
        avgSquaredQualDiff(idx) = avgSquaredQualDiff(idx) + user(i).squaredQualDiff;
        % if(user(i).starv > 0)
        %     starvCntVec(idx) = starvCntVec(idx) + 1;
        % end        
        if(user(i).starvedSegment > 0) % if starvation then value is not -1, it is greater than 0
            starvCntVec(idx) = starvCntVec(idx) + 1;
        end     
        if(user(i).visitsBadState > 0)
            vBCntVec(idx) = vBCntVec(idx) + 1;
        end
        
    end
end

if(sum(countFinishVec == 0) == 0) % a check so that there is atleast one user in each class
    avgQualitySwitchesClass = avgQualitySwitchesClass ./ countFinishVec;
    avgPrefetchTimeClass = avgPrefetchTimeClass ./ countFinishVec;
    avgDownloadTimeClass = avgDownloadTimeClass ./ countFinishVec;
    avgVideoDurationClass = avgVideoDurationClass ./ countFinishVec;
    avgQualityClass = avgQualityClass ./ countFinishVec; %in bit per sim. slot
    avgSquaredQualDiff = avgSquaredQualDiff ./ countFinishVec;
    avgQualityClass = avgQualityClass * simSlotsPerSec; % in bps
    probStarvClass = starvCntVec ./ countFinishVec;
    probVBClass = vBCntVec./countFinishVec;
    probFinishClass = countFinishVec ./ (countFinishVec + countBalkVec);
    probDropClass = countDroppedVec ./ (countFinishVec + countDroppedVec);
else
    minusOneVec = -1 * ones(1,numberOfClasses);
    avgQualitySwitchesClass = minusOneVec;
    avgPrefetchTimeClass = minusOneVec;
    avgDownloadTimeClass = minusOneVec;
    avgVideoDurationClass = minusOneVec;
    avgQualityClass = minusOneVec;
    probStarvClass = minusOneVec;
    probVBClass = minusOneVec;
    probFinishClass = minusOneVec;
    probDropClass = minusOneVec;
    avgSquaredQualDiff = minusOneVec;
end
numUsers = sum(countFinishVec);
% avgTimeinSystem = getAvgTimeinSystem(user, simSlotsPerSec);
avgTimeinSystem = -1;
% save results.mat; %For error checking
