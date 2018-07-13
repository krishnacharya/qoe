function [probStarvation, goodStates] = ...
    computeProbStarvation_zeroOrderMC_balk(tpm, badStatesIdx, probFinishing)
% compute prob. of starvation as prob. of visiting bad state in sojourn for user entering in  goodState.
% Note that prob. of starvation for user entering in badState is 1.
allStates = 1:length(badStatesIdx);
badStates = allStates(badStatesIdx == 1);
goodStates = allStates(badStatesIdx == 0);
A_inv = (eye(length(goodStates)) - tpm(goodStates, goodStates))^-1;
probStarvation = A_inv*tpm(goodStates, badStates)*probFinishing(badStates)';
% P_s = (I - P(GxG))^-1 * P(GxB) * 1(Bx1)
