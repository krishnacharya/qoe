function [WT, ExpNoOfUsers] = getWaitingTimeAnalytic(throughputVec, channelCapacity, videoQualVec,lambda, avgVidDuration, maxUsers,...
    GainVec)
    %state is number of users in the system
    %%%%%%%%%%%%%
    rho = lambda * avgVidDuration;   
    bitRateVec = zeros(1, maxUsers);
    phiVec = ones(1, maxUsers);%allocation function
    for i = 1 : maxUsers      
        [val, idx] =  min(throughputVec(i) >= videoQualVec);
        if(idx >= 2)
            bitRateVec(i) = videoQualVec(idx - 1);%choose bitrate as close to available thrpt
        else
            bitRateVec(i) = videoQualVec(1);
        end
        prodCum = 1;
        for j = 1 : i
            prodCum = prodCum * (channelCapacity / bitRateVec(j)) * GainVec(j);
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
    WT = ExpNoOfUsers / lambda;%by little's law
end













