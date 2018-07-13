function bVec = computeBufferSpacing(videoRateVec, bmin, bmax, unif)
% compute buffer thresholds for switching rate
% INPUT:
% videoRateVec: rate of available video rates
% bmin, bmax: min. and max. buffer thresholds
% unif: set to 1 if you need buffer thresholds to be uniformly spaced in [bmin, bmax]
%	0 if you need linear
% 	2 minimum to avoid non-adjacent bit rate switches
% OUTPUT: 
% bVec: vector of buffer thresholds

lmax = max(videoRateVec);
lmin = min(videoRateVec);
if (unif == 0)
    m = (bmax - bmin)/(lmax - lmin);
    c = bmax - m*lmax;
    bVec = m*videoRateVec + c;
elseif (unif == 1)
    bVec = linspace(bmin, bmax, length(videoRateVec));
elseif (unif == 2)
    bDiff = zeros(1, length(videoRateVec)-1);
    for k = 1:length(videoRateVec)-1,
        bDiff(k) = max(bDiff(k), videoRateVec(k+1)/videoRateVec(k) - 1);
        if(k < length(videoRateVec)-1)
            bDiff(k+1) = max(bDiff(k+1), 1 - videoRateVec(k)/videoRateVec(k+1));
        end
    end;
    bVec = [bmin, bmin + cumsum(bDiff)];
else
    bVec = 0;
end
end
