% INPUT parameters
arVec = 0.01; % User Arrival rate vector [lambda_1, lambda_2]
avSizeVec = 1200; % Avg. Video size 
cRate = 5e6; % Channel rate (nominal)
gammaVec = 0.0; % Vector of multipliers for balking, higher values give higher balking rate
minRateThresVec = 1e6; % Vector or minimum rate required for user satisfaction, users MAY leave if rate received is less than this
wtVec = 1; % Weight vector for max-channel rate
% A wtVec = [1, 2] and cRate of 5e6 will set max. aggregate channel rate for class 1 = 5e6 and for class 2 = 10e6.

maxUsersVec =50; % Max. number of users per class
videoRateMatrix = [0.2 0.3 0.48 0.75 1.2 1.85 2.85 4.3 5.3]*1e6; % Vector of Available bit-rates
% Row-i of videoRateMatrix gives the available bit rates for class-i.
% The vector of rates must be in ascending order.
% Append vector with zeros if the lengths are not equal

avgUsersSim = 1000; % The simulation will simulate (on average) avgUsersSim users entering the system

% DASH parameters, bmin, bmax, q_a (prefetch segments) and number of seconds per video segment respectively
bminVec = 4; bmaxVec = 10; prefetchVec = 1; secsPerSegVec = 2;
unifVec = 1; % Parameter which decides buffer thresholds, 0: linear, 1: thresholds uniformly spaced, 2: minimum required

simulationDir = '~/SCRIPTS'; % SET this to directory where you have simulation scripts (multi class with user balking)
totalSim = 3;
analysisDir = '~/SCRIPTS'; % SET this to directory where you have analysis scripts (multi class with user balking) 

cd(simulationDir);
filename = sprintf('results_50.txt');
%filename = sprintf('aR_%f_aS_%d_cR_%d_maxUsers_%d_gamma_%f.txt', ...
%    numClasses, avSizeVec, cRate, maxUsersVec, gammaVec);
fid = fopen(filename,'a+');

cd(analysisDir); % Analysis code
numClasses = length(arVec);
[pi_an, PD_an, PF_an, B_an, P_an, NQS_an, Q_an, D_an] = ...
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
printVec(fid, pi_an, length(pi_an));
fprintf(fid,'\n\n');

% START simulations here.
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
	[pi_sim, pi_2_sim, P_sim(nSim, :), B_sim(nSim, :), PF_sim(nSim, :), PD_sim(nSim, :), Q_sim(nSim, :), NQS_sim(nSim, :), D_sim(nSim, :), ...
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
	fprintf(fid,'\nSIM sanity TV (distr): %.5f', 0.5*sum(abs(pi_2_sim - pi_sim)));
	fprintf(fid,', number of users: %d', numUsers);

	fprintf(fid,'\nSIM steady state distr (SYSTEM).: ');
	printVec(fid, pi_2_sim, length(pi_2_sim));

	fprintf(fid,'\nSIM steady state distr (USER).: ');
	printVec(fid, pi_sim, length(pi_sim));
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
    fprintf(fid,'\n\n');
else
    for jj = 1:length(maxUsersVec),
        fprintf(fid,'\nSIM FINAL Res: class %d ', jj);
        fprintf(fid,'ERROR!');
        fprintf(fid,'\nSIM Dev: class %d ', jj);
        fprintf(fid,'ERROR!');
    end;
    fprintf(fid,'\n\n');
end
fclose(fid);
