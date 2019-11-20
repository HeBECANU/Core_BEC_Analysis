function GZ = g_bose(Z)
% function [GZ,err,n_ser] = g_bose(Z)
% Polylog of order 3/2 - Bose-Einstein distribution integral
% Converges for Z in [0,1) - runtime diverges as Z->1

% USER CONSTS FOR ALGORITHM CONVERGENCE
% TODO
% - optional args
% - tolerance or number of terms
% - aproximate expressions
% - consider if better to use something off file exchange
%   - https://au.mathworks.com/matlabcentral/fileexchange/23060-polylogarithm-de-jonquiere-s-function
%   - approx expansions https://github.com/wme7/Polylog
%   - https://au.mathworks.com/matlabcentral/fileexchange/37229-enhanced-computation-of-polylogarithm-aka-de-jonquieres-function
% - fix error where sumation does not terminate if any of the z vector input has not converged
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