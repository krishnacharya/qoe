function VQVec = getVideoQualVec(bVec, lmin, lmax)
    slope = (lmax - lmin) / (bVec(end) - bVec(1));
    intercept = lmax - slope * bVec(end);
    VQVec = slope * bVec + intercept;
end