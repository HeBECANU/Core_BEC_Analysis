function GZ = g_bose(Z)
% function [GZ,err,n_ser] = g_bose(Z)
% Polylog of order 3/2 - Bose-Einstein distribution integral
% Converges for Z in [0,1) - runtime diverges as Z->1

% USER CONSTS FOR ALGORITHM CONVERGENCE
tol_err=1e-20;   % incremental error from evaluating series sum

zmask = Z >= 0.9999;
if any(zmask)
    warning('g_bose will not converge for Z>=1 !')
    Z(zmask) = nan;
end
% else

GZ=0;   % initialise output
err=inf;    % initialise incremental error

n_ser=0;    % summand order
while err>tol_err
    n_ser=n_ser+1;
    GZ_temp=GZ;     % save last order series sum
    
    % evaluate next series term
    GZ=GZ+(Z.^n_ser)/(n_ser^(3/2));

    % evaluate error from this term
    if n_ser>1  % skip evaluateion (div) for first iter
        err=max(abs(GZ-GZ_temp)./GZ_temp);     % incremental fractional diff
    end
end
end