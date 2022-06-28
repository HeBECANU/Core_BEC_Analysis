
function gauss_amp=gauss_2d_rot_function(cord, sigma, center, theta, offset, amp,convention)
% gaussian function for fitting beam profiles
%
% Known BUGS/ Possible Improvements
%   - add optional input 'normalization',['peak','sum','none']
%
% Author: Bryce Henson
% email: Bryce.m.henson+gaussian_function_2d@gmail.com
% Last revision:2019-03-07

%------------- BEGIN CODE --------------
if isempty(center)
    center=[0,0];
end
if isempty(theta)
    theta=0;
end

if isempty(convention)
    convention='sigma';
end

if isempty(amp)
    amp=1;
end

% this function folows the convention of 
%https://en.wikipedia.org/wiki/Gaussian_beam
% with exp(-2 r^2 / w^2 )

xvals=cord(:,1);
yvals=cord(:,2);
xc      = center(1);
yc      = center(2);
theta   = (theta/180)*pi;
xm      = (xvals-xc).*cos(theta) - (yvals-yc).*sin(theta);
ym      = (yvals-xc).*sin(theta) + (yvals-yc).*cos(theta);

if strcmp(convention,'sigma')
    gauss_amp   = exp(-(1/2)*((xm./sigma(1)).^2 + (ym./sigma(2)).^2));
elseif strcmp(convention,'waist')
    gauss_amp   = exp(-(2)*((xm./sigma(1)).^2 + (ym./sigma(2)).^2));
end


gauss_amp=gauss_amp*amp;

gauss_amp=gauss_amp+offset;

end
