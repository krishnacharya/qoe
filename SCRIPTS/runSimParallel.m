%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 0: Set INPUT parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Arrival rate vector [lambda_1, lambda_2, ...] where \lambda_k = class-k 
% arrival rate.
% arVec = 0.01; % E.g. [0.01, 0.01]
lambdas = 0.01 : 0.01: 0.2;% arrival rates
nOfLambdas = length(lambdas);
% Avg. Video size [\theta_1, \theta_2, ...]
avSizeVec = 200; % E.g. [1200, 600]
% Channel rate (nominal) [aggregate rate available for class K  = wtVec[k]*cRate]
%cRate = 5e6; % E.g., 5e6, 5Mbps

% Vector of multipliers for balking, higher values give higher balking rate
gammaVec = 0; % E.g., [100.0, 100.0]

% Vector or minimum rate required for user satisfaction, users MAY leave if 
% rate received is less than this. [R_1, R_2, ...] R_k = rate threshold for 
% class-k user, people are not satisfied if they receive video-bit rate less 
% than R_k.  
minRateThresVec = 1e6; % E.g., [5e5, 1e6]

% Weight vector for max-channel rate
% A wtVec = [1, 2] and cRate of 5e6 will set max. aggregate channel rate for 
% class 1 = 5e6 and for class 2 = 10e6.
%wtVec = 1; % E.g., [1, 2]

% Max. number of users per class [N_1, N_2, ...]
% N_k = max. number of users of class k that can be there in the system.
maxUsersVec = 20; % E.g., [10, 20]

% Matrix of Available bit-rates
% Row-i of videoRateMatrix gives the available bit rates for class-i.
% The vector of rates must be in ascending order.
% Append vector with zeros if the lengths are not equal
videoRateMatrix = [0.2, 0.3, 0.48, 0.75, 1.2, 1.85, 2.85, 4.3, 5.3] * 1e6; 

% The simulation will simulate (on average) avgUsersSim users entering the system
avgUsersSim = 1000;% E.g. 2000

% DASH parameters, bmin, bmax, q_a (prefetch segments) and number of seconds per 
% video segment respectively
bminVec = 4; bmaxVec = 10; prefetchVec = 2; secsPerSegVec = 2; % Vector valued

% Parameter which decides buffer thresholds, 0: linear, 1: thresholds uniformly 
% spaced, 2: minimum required
unifVec = 1; % vector valued

% For throughput calculations
channelStatesVec = [2, 3, 5] * 1e6;% 1Mbps, 3 Mbps and 5 Mbps
channelStatesDistr = [0.2, 0.5, 0.3];%must sum to 1, prob distr
throughputVec = zeros(1, maxUsersVec); % currently for 1 class
channelCapacity = channelStatesVec * channelStatesDistr'; %channel state when only 1 user can run R
alpha=1; beta=1; %exponenets used in priority calculation
for i = 1 : maxUsersVec
    throughputVec(i) = mean(getThroughput(i, channelStatesVec, channelStatesDistr,alpha, beta));%throughput when there are i users in the system
end
[Glimit, GainVec] = getGainLimit(channelStatesVec, channelStatesDistr, alpha, beta, maxUsersVec);%independent of lamda
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
WaitingTimeVecAnal = zeros(1, nOfLambdas);
ENofUVecAnal = zeros(1, nOfLambdas);
WaitingTimeVecSim = zeros(1, nOfLambdas);
ENofUVecSim = zeros(1, nOfLambdas);
StabilityRatioVec = zeros(1,nOfLambdas);

piSimVec = cell(1, nOfLambdas); pi_2SimVec = cell(1, nOfLambdas); probStarvSimVec = zeros(1, nOfLambdas); probDropClassSimVec = zeros(1, nOfLambdas);
probVBSimVec = zeros(1, nOfLambdas); probFinishSimVec = zeros(1, nOfLambdas); prefetchDelaySimVec = zeros(1, nOfLambdas);
avgQualitySwitchesSimVec = zeros(1, nOfLambdas); avgQualitySimVec = zeros(1, nOfLambdas), PrefetchDelayijSimVec = cell(1, nOfLambdas);

piAnalVec = cell(1, nOfLambdas); pyeDirectAnalVec = cell(1, nOfLambdas); probBlockingAnalVec = zeros(1, nOfLambdas); probFinishAnalVec = zeros(1, nOfLambdas);
probVBAnalVec = zeros(1, nOfLambdas); probStarvAnalVec = zeros(1, nOfLambdas); prefetchDelayAnalVec = zeros(1, nOfLambdas);
avgQualitySwitchesAnalVec = zeros(1, nOfLambdas); avgQualityAnalVec = zeros(1, nOfLambdas), PrefetchDelayijAnalVec = cell(1, nOfLambdas);

parfor i = 1:nOfLambdas
    %Simulation
    %M('lambda') = lambdas(i);
    [piSimVec{i}, pi_2SimVec{i}, userListSim, avgTimeinSystem, probStarvSimVec(i), probVBSimVec(i), probFinishSimVec(i), probDropClassSimVec(i), ...
    avgQualitySimVec(i), avgQualitySwitchesSimVec(i),prefetchDelaySimVec(i),PrefetchDelayijSimVec{i}, FreqMatrix, avgDownloadTimeSim, avgVideoDurationSim, numUsers, simTime] = ... 
    simscript2(lambdas(i), prefetchVec, avSizeVec, secsPerSegVec, gammaVec, minRateThresVec, ...
    maxUsersVec, videoRateMatrix, unifVec, bminVec, bmaxVec, avgUsersSim, throughputVec) 
    ExN = 0;
    for j = 1 : 21
        ExN = ExN + (j-1) * piSimVec{i}(j);
    end
    ENofUVecSim(i) = ExN;
    b = videoRateMatrix;%used to get minimum video quality
    StabilityRatioVec(i) = lambdas(i) * avSizeVec * b(1) / (channelCapacity * Glimit);
    %Analytic part
    [piAnalVec{i}, pyeDirectAnalVec{i}, probBlockingAnalVec(i), probFinishAnalVec(i), probVBAnalVec(i), probStarvAnalVec(i), avgQualitySwitchesAnalVec(i), ...
    avgQualityAnalVec(i), prefetchDelayAnalVec(i), PrefetchDelayijAnalVec{i}] = analyticUserMC(lambdas(i), avSizeVec, gammaVec, minRateThresVec, maxUsersVec,...
    videoRateMatrix, prefetchVec, secsPerSegVec, bminVec, bmaxVec, unifVec, throughputVec) 
    %TO CHANGE%%%%
    [WaitingTimeVecAnal(i), ENofUVecAnal(i)] = getWaitingTimeAnalytic(throughputVec, channelCapacity, videoRateMatrix, lambdas(i), avSizeVec, maxUsersVec,...
    GainVec);
end
WaitingTimeVecSim = ENofUVecSim ./ lambdas;

% fprintf(fid,'Pi user perspective\n');
% printVec(fid, pi, length(pi));

% fprintf(fid,'Pi per slot glanurity\n');
% printVec(fid, pi_2, length(pi_2));
