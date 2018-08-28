function [WT, ExpNoOfUsers] = getWaitingTimeAnalytic(throughputVec, channelCapacity, videoQualVec, lamda, avgVidDuration, maxUsers)
    %state is number of users in the system
    rho = lamda * avgVidDuration;
    thrptVec = zeros(1, maxUsers);%throughput vector
    bitRateVec = zeros(1, maxUsers);
    phiVec = ones(1, maxUsers);%allocation function
    for i = 1 : maxUsers      
        [val, idx] =  min(throughputVec(i) > videoQualVec);
        bitRateVec(i) = videoQualVec(idx - 1);%choose bitrate as close to available thrpt
        prodCum = 1;
        for j = 1 : i
            prodCum = prodCum * (channelCapacity / bitRateVec(j));
        end
        phiVec(i) = prodCum;
    end
    Jm = 0; % normalization term
    for i = 1 : maxUsers
        Jm = Jm + rho^i / phiVec(i);
    end
    ExpNoOfUsers = 0;%expected number of users in system
    %Prob(N=n) = (1/Jm) * rho^n / phiVec(n) 
    for i = 1: maxUsers
        ExpNoOfUsers = ExpNoOfUsers + i * rho^i / phiVec(i);
    end
    ExpNoOfUsers = ExpNoOfUsers / Jm;
    WT = ExpNoOfUsers / lamda;%by little's law
end













