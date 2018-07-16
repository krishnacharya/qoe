%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 0: Set INPUT parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Arrival rate vector [lambda_1, lambda_2, ...] where \lambda_k = class-k 
% arrival rate.
arVec = 0.01; % E.g. [0.01, 0.01]

% Avg. Video size [\theta_1, \theta_2, ...]
avSizeVec = 1200; % E.g. [1200, 600
% Channel rate (nominal) [aggregate rate available for class K  = wtVec[k]*cRate]
cRate = 5e6; % E.g., 5e6, 5Mbps

% Vector of multipliers for balking, higher values give higher balking rate
gammaVec = 100.0; % E.g., [100.0, 100.0]

% Vector or minimum rate required for user satisfaction, users MAY leave if 
% rate received is less than this. [R_1, R_2, ...] R_k = rate threshold for 
% class-k user, people are not satisfied if they receive video-bit rate less 
% than R_k.  
minRateThresVec = 1e6; % E.g., [5e5, 1e6]

% Weight vector for max-channel rate
% A wtVec = [1, 2] and cRate of 5e6 will set max. aggregate channel rate for 
% class 1 = 5e6 and for class 2 = 10e6.
wtVec = 1; % E.g., [1, 2]

% Max. number of users per class [N_1, N_2, ...]
% N_k = max. number of users of class k that can be there in the system.
maxUsersVec = 100; % E.g., [10, 20]

% Matrix of Available bit-rates
% Row-i of videoRateMatrix gives the available bit rates for class-i.
% The vector of rates must be in ascending order.
% Append vector with zeros if the lengths are not equal
videoRateMatrix = [0.2, 0.3, 0.48, 0.75, 1.2, 1.85, 2.85, 4.3, 5.3] * 1e6; 

% The simulation will simulate (on average) avgUsersSim users entering the system
avgUsersSim = 1000; % E.g. 2000

% DASH parameters, bmin, bmax, q_a (prefetch segments) and number of seconds per 
% video segment respectively
bminVec = 4; bmaxVec = 10; prefetchVec = 1; secsPerSegVec = 2; % Vector valued

% Parameter which decides buffer thresholds, 0: linear, 1: thresholds uniformly 
% spaced, 2: minimum required
unifVec = 1; % vector valued

% SET this to directory where you have simulation scripts (multi class with user balking)
simulationDir = '~/SCRIPTS';

% Total number of simulations that you wish to run.
totalSim = 3;

% SET this to directory where you have analysis scripts (multi class with user balking) 
analysisDir = '~/SCRIPTS';

filename = 'firstRun.txt';
fid = fopen(filename,'a+'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pi, pi_2, user, probStarvClass, probVBClass, probFinishClass, probDropClass, avgQualityClass, avgQualitySwitchesClass, ...
    avgPrefetchTimeClass, avgDownloadTimeClass, avgVideoDurationClass, numUsers, simTime] = ...
    simscript(arVec, avSizeVec, cRate, gammaVec, minRateThresVec, wtVec, maxUsersVec,videoRateMatrix, unifVec, avgUsersSim);

[pye, probBlocking, probFinishing, probVB, probStarvation, avgQualSwitches, ...
    avgQuality, prefetchDelay] = ...
    userMC_firstOrderMC_PF_balk(arVec, avSizeVec, ...
    cRate, gammaVec, minRateThresVec, wtVec, maxUsersVec, ...
    videoRateMatrix, prefetchVec, secsPerSegVec, bminVec, bmaxVec, unifVec);

% fprintf(fid,'Pi user perspective\n');
% printVec(fid, pi, length(pi));

% fprintf(fid,'Pi per slot glanurity\n');
% printVec(fid, pi_2, length(pi_2));
