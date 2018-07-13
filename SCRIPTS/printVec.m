function printVec(fid, vector, length)
% PRINT first 'length' entries of 'vector' into file with file descriptor 'fid'.
for ii=1:length,
    if(abs(vector(ii)) <= 1e-5)
        fprintf(fid, '0, ');
    elseif (abs(vector(ii) - 1) <= 1e-5)
        fprintf(fid, '1, ');
    else
        fprintf(fid, '%.5f, ',  vector(ii));
    end;
end;
