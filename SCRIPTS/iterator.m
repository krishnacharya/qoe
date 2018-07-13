function [V, D] = iterator(TPM)
% Find stationary distribution by using iteration, i.e., use
% pi[k] P^k --> pi^* where pi^* is the stationary distribution.
pi(1:length(TPM)) = 1/length(TPM);
error = 1;
budget = 1000;
while(error >= 1e-10)
    pi2 = pi*TPM;
    error = sum((pi - pi2).^2);
    pi = pi2;
    budget = budget - 1;
    if(budget == 0),
        break;
    end;    
end;
V = pi';
D = 1;
