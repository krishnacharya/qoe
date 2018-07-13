function tpm = trm2Embeddedtpm(trm)
% Get the embedded TPM corresponding to rate matrix trm.
tpm = zeros(size(trm));
for ii = 1:length(trm(:,1))
    if(abs(trm(ii,ii)) > 0)
        tpm(ii,:) = abs(trm(ii,:)/trm(ii,ii));
        tpm(ii,ii) = 0;
    else
        tpm(ii,ii) = 1;
    end
end
