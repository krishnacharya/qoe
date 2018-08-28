function GLimit = getGainLimit(M)
% GLimit or G*, is required for verfying the stability of the system for
% certain values of lamda, pho < R G* / bmin
R = M('channelCapacity');
prev = 0;
i = 1;
nplot = zeros(1, 100);
    while true
        M('n') = i;
        nplot(i) = i * mean(getThroughput(M)) / R;
        if(abs(nplot(i) - prev) / prev <= 0.001) % stop if less than 0.1% change
            break;
        end
        prev = nplot(i);
%         disp(i);
%         disp(nplot(i));
        i = i+1;
    end
GLimit = prev;
end