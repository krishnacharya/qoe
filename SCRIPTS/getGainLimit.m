function [GLimit,GainVec] = getGainLimit(M)
% GLimit or G*, is required for verfying the stability of the system for
% certain values of lamda, pho < R G* / bmin
R = M('channelCapacity');
maxUserVec = M('maxUsersVec');
prev = 0;
i = 1;
GainVec = zeros(1, maxUserVec);
    while true
        M('n') = i;
        GainVec(i) = i * mean(getThroughput(M)) / R;
        if((abs(GainVec(i) - prev) / prev <= 0.001) && (i > maxUserVec)) % stop if less than 0.1% change
            break;
        end
        prev = GainVec(i);
%         disp(i);
%         disp(GainVec(i));
        i = i+1;
    end
GLimit = prev;
end