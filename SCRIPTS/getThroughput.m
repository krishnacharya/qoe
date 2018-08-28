function usrAvgThrptVec = getThroughput(M)
%get the throughput for a user when there are n users in the system.
%for now assume we have only 1 class, so each of their average throughputs
%should be the same in the long run

% channelStatesVec = [2, 3, 5] * 1e6;% 1Mbps, 3 Mbps and 5 Mbps
% channelStatesDistr = [0.2, 0.5, 0.3];%must sum to 1, prob distr

% picking relevant variable from input dictionary
n = M('n');%number of users e.g PF(10)
channelStatesVec = M('channelStatesVec');
channelStatesDistr = M('channelStatesDistr');
alpha = M('alpha');
beta = M('beta');

slots = 1e5;
usrAvgThrptVec = ones(1, n) * 1e6;%Avg Throughput vector, update after each slot
% alpha=1;
% beta=1;%for prop sharing in LTE
    for i = 1 : slots
        usrInstThrptVec = zeros(1, n);
        for j = 1 : n
           channelType = selectUserType(channelStatesDistr);
           usrInstThrptVec(j) = channelStatesVec(channelType);
        end
        priorityVec = (usrInstThrptVec .^ alpha) ./ (usrAvgThrptVec .^ beta);
        maxIndices = find(priorityVec == max(priorityVec));%all the tied max Priorities
        chosenIdx = maxIndices(randsample(length(maxIndices), 1));
        temp = usrAvgThrptVec(chosenIdx);
        usrAvgThrptVec = usrAvgThrptVec * (i-1) / i;
        usrAvgThrptVec(chosenIdx) = (temp * (i-1) + usrInstThrptVec(chosenIdx)) / i;
    end
end