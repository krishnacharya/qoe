%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 0: Set INPUT parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Arrival rate vector [lambda_1, lambda_2, ...] where \lambda_k = class-k 
% arrival rate.
% arVec = 0.01; % E.g. [0.01, 0.01]
M = containers.Map();
lambdas = 0.01 : 0.05: 0.2;% arrival rates
nOfLamdas = length(lambdas);
M('lamda') = 0.01;
% Avg. Video size [\theta_1, \theta_2, ...]
M('avSizeVec') = 1200; % E.g. [1200, 600]
% Channel rate (nominal) [aggregate rate available for class K  = wtVec[k]*cRate]
%cRate = 5e6; % E.g., 5e6, 5Mbps

% Vector of multipliers for balking, higher values give higher balking rate
M('gammaVec') = 0; % E.g., [100.0, 100.0]

% Vector or minimum rate required for user satisfaction, users MAY leave if 
% rate received is less than this. [R_1, R_2, ...] R_k = rate threshold for 
% class-k user, people are not satisfied if they receive video-bit rate less 
% than R_k.  
M('minRateThresVec') = 1e6; % E.g., [5e5, 1e6]

% Weight vector for max-channel rate
% A wtVec = [1, 2] and cRate of 5e6 will set max. aggregate channel rate for 
% class 1 = 5e6 and for class 2 = 10e6.
%wtVec = 1; % E.g., [1, 2]

% Max. number of users per class [N_1, N_2, ...]
% N_k = max. number of users of class k that can be there in the system.
M('maxUsersVec') = 20; % E.g., [10, 20]

% Matrix of Available bit-rates
% Row-i of videoRateMatrix gives the available bit rates for class-i.
% The vector of rates must be in ascending order.
% Append vector with zeros if the lengths are not equal
M('videoRateMatrix') = [0.2, 0.3, 0.48, 0.75, 1.2, 1.85, 2.85, 4.3, 5.3] * 1e6; 

% The simulation will simulate (on average) avgUsersSim users entering the system
M('avgUsersSim') = 100; % E.g. 2000

% DASH parameters, bmin, bmax, q_a (prefetch segments) and number of seconds per 
% video segment respectively
M('bminVec') = 4; M('bmaxVec') = 10; M('prefetchVec') = 2; M('secsPerSegVec') = 2; % Vector valued

% Parameter which decides buffer thresholds, 0: linear, 1: thresholds uniformly 
% spaced, 2: minimum required
M('unifVec') = 1; % vector valued

% For throughput calculations
M('channelStatesVec') = [2, 3, 5] * 1e6;% 1Mbps, 3 Mbps and 5 Mbps
M('channelStatesDistr') = [0.2, 0.5, 0.3];%must sum to 1, prob distr
throughputVec = zeros(1, M('maxUsersVec')); % currently for 1 class
M('channelCapacity') = channelStatesVec * channelStatesDistr'; %channel state when only 1 user can run R
M('alpha')=1; M('beta')=1; %exponenets used in priority calculation
for i = 1 : M('maxUsersVec')
    M('n') = i;
    throughputVec(i) = mean(getThroughput(M));%throughput when there are i users in the system
end
M('Glimit') = getGainLimit(M);%independent of lamda
M('throughputVec') = throughputVec;
% % SET this to directory where you have simulation scripts (multi class with user balking)
% simulationDir = '~/SCRIPTS';
% 
% % Total number of simulations that you wish to run.
% totalSim = 3;
% 
% % SET this to directory where you have analysis scripts (multi class with user balking) 
% analysisDir = '~/SCRIPTS';
% 
% filename = 'firstRun.txt';
% fid = fopen(filename,'a+'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1 : nOfLamdas
%     M('lamda') = lambdas(i);
%     [piSim, pi_2Sim, userListSim, avgTimeinSystem, probStarvSim, probVBSim, probFinishSim, probDropClass, avgQualitySim, avgQualitySwitchesSim,prefetchDelaySim,PrefetchDelayijSim, FreqMatrix, avgDownloadTimeSim, avgVideoDurationSim, numUsers, simTime] = ... 
%     simscript(M);
% 
%     [piAnal, pyeDirectAnal, probBlockingAnal, probFinishAnal, probVBAnal, probStarAnal, avgQualSwitchesAnal, ...
%         avgQualityAnal, prefetchDelayAnal, prefetchDelayijAnal] =   analyticUserMC(M);
% end
WTVecAnal = zeros(1, nOfLamdas);
ENofUVecAnal = zeros(1, nOfLamdas);
WTVecSim = zeros(1, nOfLamdas);
ENofUVecSim = zeros(1, nOfLamdas);
StabilityRatio = zeros(1,nOfLamdas);
for i = 1:nOfLamdas
    %Simulation
    M('lamda') = lambdas(i);
    [piSim, pi_2Sim, userListSim, avgTimeinSystem, probStarvSim, probVBSim, probFinishSim, probDropClass, avgQualitySim, avgQualitySwitchesSim,prefetchDelaySim,PrefetchDelayijSim, FreqMatrix, avgDownloadTimeSim, avgVideoDurationSim, numUsers, simTime] = ... 
    simscript2(M);   
    ExN = 0;
    for j = 1 : 21
        ExN = ExN + (j-1) * pi_2Sim(j);
    end
    ENofUVecSim(i) = ExN;
    b = M('videoRateMatrix');%used to get minimum video quality
    StabilityRatio(i) = lambdas(i) * M('avSizeVec') * b(1) / (M('channelCapacity') * M('Glimit'));
    %Analytic part
    [piAnal, pyeDirectAnal, probBlockingAnal, probFinishAnal, probVBAnal, probStarAnal, avgQualSwitchesAnal, ...
    avgQualityAnal, prefetchDelayAnal, prefetchDelayijAnal] =   analyticUserMC(M);   
    %TO CHANGE%%%%
    %[WTVecAnal(i), ENofUVecAnal(i)] = getWaitingTimeAnalytic(throughputVec, channelRate, videoQualVec, lambdas(i), avgVidDuration, maxUsers);
end
WTVecSim = ENofUVecSim ./ lambdas;



% fprintf(fid,'Pi user perspective\n');
% printVec(fid, pi, length(pi));

% fprintf(fid,'Pi per slot glanurity\n');
% printVec(fid, pi_2, length(pi_2));
