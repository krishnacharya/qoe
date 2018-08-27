function GLimit = getGainLimit(R)
% GLimit or G*, is required for verfying the stability of the system for
% certain values of lamda, pho < R G* / bmin
prev = 0;
i = 1;
    while true
        nplot(i) = i * mean(getThroughput(i)) / R;
        if(abs(nplot(i) - prev) / prev <= 0.001)
            break;
        end
        prev = nplot(i);
        disp(i);
        disp(nplot(i));
        i = i+1;
    end
GLimit = prev;
end