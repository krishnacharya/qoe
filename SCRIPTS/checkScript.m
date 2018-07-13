%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 0: Set INPUT parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Arrival rate vector [lambda_1, lambda_2, ...] where \lambda_k = class-k 
% arrival rate.
arVec = 0.01; % E.g. [0.01, 0.01]

% Avg. Video size [\theta_1, \theta_2, ...]
avSizeVec = 1200; % E.g. [1200, 600]

% Channel rate (nominal) [aggregate rate available for class K  = wtVec[k]*cRate]
cRate = 5e6; % E.g., 5e6

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
maxUsersVec = 20; % E.g., [10, 20]

% Matrix of Available bit-rates
% Row-i of videoRateMatrix gives the available bit rates for class-i.
% The vector of rates must be in ascending order.
% Append vector with zeros if the lengths are not equal
videoRateMatrix = [0.2 0.3 0.48 0.75 1.2 1.85 2.85 4.3 5.3]*1e6; 
% E.g., [0.2 0.3 4.3 5.3; 0.2 0.3 2.3 5.3]*1e6

% The simulation will simulate (on average) avgUsersSim users entering the system
avgUsersSim = 2000; % E.g. 2000


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: RUN Analysis code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(simulationDir);
sim_time = clock;
filename = sprintf('results_%02d-%02d-%02d-%02d-%02d.txt', sim_time(1:5));
fid = fopen(filename,'a+');

% Print all the simulation parameters into results file.
fprintf(fid,'\n\nRESULTS for ');
for ii = 1:length(arVec),
    fprintf(fid,'Class %d: aR: %.3f, vSize: %d, ', ii, arVec(ii), avSizeVec(ii));
    fprintf(fid,'gamma: %.1f, min. rate thres.: %.1f, ', gammaVec(ii), ...
        minRateThresVec(ii));
    fprintf(fid,'weight: %d, max. users: %d, ', wtVec(ii), maxUsersVec(ii));
    fprintf(fid,'bmin: %.1f, bmax: %.1f, ', bminVec(ii), bmaxVec(ii));
    fprintf(fid,'prefetch thres.: %.1f, secs/seg: %.1f, ', prefetchVec(ii), ...
        secsPerSegVec(ii));
    fprintf(fid,'unif: %d\n', unifVec(ii));
end;
fprintf(fid,'cRate: %.1f\n', cRate);




cd(analysisDir); % Analysis code in analysisDir
numClasses = length(arVec);
[distr_an, PD_an, PF_an, B_an, P_an, NQS_an, Q_an, D_an] = ...
    userMC_firstOrderMC_PF_balk(arVec, avSizeVec, cRate, gammaVec, ...
    minRateThresVec, wtVec, maxUsersVec, videoRateMatrix, prefetchVec, ...
    secsPerSegVec, bminVec, bmaxVec, unifVec);

% Print analysis results
cd(simulationDir);
RQS_an =  NQS_an ./ avSizeVec;
for jj = 1:length(maxUsersVec),    
    fprintf(fid,'\nAnalysis Res: class %d ', jj);
    fprintf(fid,'S: %.5f, B: %.5f, Q: %.3f, NQS: %.3f, RQS: %.3f, D: %.3f, PD: %.3f, PF: %.5f,', ...
        P_an(jj), B_an(jj), Q_an(jj), NQS_an(jj), ...
        RQS_an(jj), D_an(jj), PD_an(jj), PF_an(jj));
    
end;
fprintf(fid,'\nAnalysis steady state distr (USER).: ');
printVec(fid, distr_an, length(distr_an));
fprintf(fid,'\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2 : RUN Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(simulationDir);
P_sim = zeros(totalSim,numClasses);
B_sim = zeros(totalSim,numClasses);
Q_sim = zeros(totalSim,numClasses);
NQS_sim = zeros(totalSim,numClasses);
D_sim = zeros(totalSim,numClasses);
DT_sim = zeros(totalSim,numClasses);
VD_sim = zeros(totalSim,numClasses);
RQS_sim = zeros(totalSim,numClasses);
PD_sim = zeros(totalSim,numClasses);

for nSim = 1:totalSim,
    sim_start = clock;
    fprintf(fid,'\nSIM %d START: [%02d-%02d-%02d], %02d:%02d', nSim, ...
        sim_start(1:5));

	[distr_sim, distr_2_sim, P_sim(nSim, :), B_sim(nSim, :), PF_sim(nSim, :), ...
	    PD_sim(nSim, :), Q_sim(nSim, :), NQS_sim(nSim, :), D_sim(nSim, :), ...
    	DT_sim(nSim, :), VD_sim(nSim, :), numUsers, simTime] = ...
	    	script_Imen_SP_MC_PF_balk(arVec, avSizeVec, ...
		        cRate, gammaVec, minRateThresVec, wtVec, maxUsersVec, ...
    		    videoRateMatrix, unifVec, avgUsersSim);
	% Print simulation results
	RQS_sim =  NQS_sim ./ VD_sim;
	for jj = 1:length(maxUsersVec),
	    fprintf(fid,'\nSIM Res: class %d ', jj);
	    fprintf(fid,'S: %.5f, B: %.5f, Q: %.3f, NQS: %.3f, RQS: %.3f, D: %.3f, PD: %.3f, PF: %.5f, DT: %.3f', ...
        	P_sim(nSim, jj), B_sim(nSim, jj), Q_sim(nSim,jj), NQS_sim(nSim,jj), ...
	        RQS_sim(nSim,jj), D_sim(nSim,jj), PD_sim(nSim,jj), PF_sim(nSim,jj), DT_sim(nSim,jj));
    
	end;
	fprintf(fid,'\nSIM duration: %.1f', simTime);
	fprintf(fid,'\nSIM sanity TV (distr): %.5f', 0.5*sum(abs(distr_2_sim - distr_sim)));
	fprintf(fid,', number of users: %d', numUsers);

	fprintf(fid,'\nSIM steady state distr (SYSTEM).: ');
	printVec(fid, distr_2_sim, length(distr_2_sim));

	fprintf(fid,'\nSIM steady state distr (USER).: ');
	printVec(fid, distr_sim, length(distr_sim));
	
    sim_stop = clock;
    fprintf(fid,'\nSIM %d END: [%02d-%02d-%02d], %02d:%02d\n\n', nSim, ...
        sim_stop(1:5));
	fprintf(fid,'\n\n');
end;

%% PRINT FINAL Simulation results
if (sum(P_sim < 0) == 0)
    for jj = 1:length(maxUsersVec),
        fprintf(fid,'\nSIM FINAL Res: class %d ', jj);
        fprintf(fid,'S: %.5f, B: %.5f, Q: %.3f, NQS: %.3f, RQS: %.3f, D: %.3f, PD: %.3f', ...
            mean(P_sim(:, jj)), mean(B_sim(:, jj)), mean(Q_sim(:,jj)), ...
            mean(NQS_sim(:,jj)), mean(RQS_sim(:,jj)), mean(D_sim(:,jj)), ...
            mean(PD_sim(:,jj)));
        fprintf(fid,'\nSIM Dev: class %d ', jj);
        fprintf(fid,'S: %.5f, B: %.5f, Q: %.3f, NQS: %.3f, RQS: %.3f, D: %.3f, PD: %.3f', ...
            std(P_sim(:, jj)), std(B_sim(:, jj)), std(Q_sim(:,jj)), ...
            std(NQS_sim(:,jj)), std(RQS_sim(:,jj)), std(D_sim(:,jj)), ...
            std(PD_sim(:,jj)));
    end;
else
    for jj = 1:length(maxUsersVec),
        fprintf(fid,'\nSIM FINAL Res: class %d ', jj);
        fprintf(fid,'ERROR!');
        fprintf(fid,'\nSIM Dev: class %d ', jj);
        fprintf(fid,'ERROR!');
    end;
end
sim_done = clock;
fprintf(fid,'\nSIM DONE: [%02d-%02d-%02d], %02d:%02d\n\n', ...
        sim_done(1:5));
fprintf(fid,'\n\n');
fclose(fid);
