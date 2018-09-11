function [GLimit,GainVec] = getGainLimit(channelStatesVec, channelStatesDistr,alpha, beta, maxUsersVec, videoRateMatrix)
% GLimit or G*, is required for verfying the stability of the system for
% certain values of lamda, pho < R G* / bmin
R = channelStatesDistr * channelStatesVec';
prev = 0;
i = 1;
GainVec = zeros(1, maxUsersVec);
    while true
        GainVec(i) = i * getThroughput(i, channelStatesVec, channelStatesDistr, alpha, beta, videoRateMatrix) / R;
        if((abs(GainVec(i) - prev) / prev <= 0.001) && (i > maxUsersVec)) % stop if less than 0.1% change
            break;
        end
        prev = GainVec(i);
        disp(i);
        disp(GainVec(i));
%         disp(GainVec(i));
        i = i + 1;
    end
GLimit = prev;
end