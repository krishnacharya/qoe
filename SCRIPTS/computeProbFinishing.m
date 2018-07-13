function probFinishing = computeProbFinishing(trm2, badStates)
% compute prob. of user actually finishing service without leaving system.
listStates = 1:length(badStates);% this is  just all the states in the markov chain
validStates = listStates(badStates == 0 | badStates == 1);
finishState = listStates(badStates == 2);
emTpm = trm2Embeddedtpm(trm2);
P = emTpm(validStates, validStates);
I = eye(length(validStates));
probFinishing = (I - P)^-1*(emTpm(validStates, finishState));
probFinishing = probFinishing';
