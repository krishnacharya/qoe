lda = 0.01:0.02:0.2;
n = length(lda);
ExNVec = zeros(1, n);

avSizeVec = 1200; % E.g. [1200, 600]
cRate = 5e6; % E.g., 5e6, 5Mbps
gammaVec = 0; % E.g., [100.0, 100.0]
minRateThresVec = 1e6; % E.g., [5e5, 1e6]
wtVec = 1; % E.g., [1, 2]
maxUsersVec = 20; % E.g., [10, 20]
videoRateMatrix = [0.2, 0.3, 0.48, 0.75, 1.2, 1.85, 2.85, 4.3, 5.3] * 1e6; 
avgUsersSim = 100; % E.g. 2000
bminVec = 4; bmaxVec = 10; prefetchVec = 2; secsPerSegVec = 2; % Vector valued
unifVec = 1; % vector valued

for i = 1:n
    [piSim, pi_2Sim, userListSim, avgTimeinSystem, probStarvSim, probVBSim, probFinishSim, probDropClass, avgQualitySim, avgQualitySwitchesSim,prefetchDelaySim,PrefetchDelayijSim, FreqMatrix, avgDownloadTimeSim, avgVideoDurationSim, numUsers, simTime] = ... 
simscript(lda(i), prefetchVec, avSizeVec, secsPerSegVec, cRate, gammaVec, minRateThresVec, wtVec, maxUsersVec, videoRateMatrix, unifVec, bminVec, bmaxVec, avgUsersSim);
    ExN = 0;
    for j = 1 : 21
        ExN = ExN + (j-1)*pi_2Sim(j)
    end
    ExNVec(i) = ExN;
end
