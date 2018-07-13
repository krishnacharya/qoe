function tpm = trm2tpm(trm)
% Get TPM (first-order approx./ uniformized DTMC) corr. to rate matrix trm
D = abs(diag(trm));
h = 0.1*min(1./(D(D > 0)));
tpm = eye(length(trm(1,:))) + trm*h;
