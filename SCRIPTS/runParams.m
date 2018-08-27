channelRate = 5e6;
videoQualVec = [0.2, 0.3, 0.48, 0.75, 1.2, 1.85, 2.85, 4.3, 5.3] * 1e6;
avgVidDuration = 1200;
maxUsers = 20;
lamdaVec = 0.01 : 0.02: 0.2;
n = length(lamdaVec);
WT = zeros(1, n);
ENofU = zeros(1, n);
for i = 1 : n
    [WT(i), ENofU(i)] = getWaitingTimeAnalytic(channelRate, videoQualVec, lamdaVec(i), avgVidDuration, maxUsers);
end
