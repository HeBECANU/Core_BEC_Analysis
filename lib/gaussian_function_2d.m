
function gauss=gaussian_function_2d(filter_size, sigma)
% gaussian_function_2d - Compute a two dimensional gaussian with different widths in x,y
% Matlabs inbilt gaussian function fspecial('gaussian', [30 30], 4) has a the same sd in each dirn
% which makes it inapropriate to use as a blur kernel for 'squished' 2d plots.
% does not inculde rotation
%
% Syntax:  [data,import_opts]=import_data(import_opts)
%
% Inputs:
%    filter_size	- vector[1x2], of the output matrix
%    sigma          - vector[1x2], width of the gaussian in units of the filter_size
%
% Outputs:
%    gauss      - matrix, sampled gaussian
%
%
% Example: 
%   imagesc(Gaussian_filter2d([100 100],[10 10]))
%   set(gca,'dataAspectRatio',[1 1 1])

% Other m-files required: none
% Also See: none
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%   - unit testing
%
% Author: Bryce Henson
% email: Bryce.m.henson+gaussian_function_2d@gmail.com
% Last revision:2019-03-07

%------------- BEGIN CODE --------------

sizex=filter_size(1);
sizey=filter_size(1);
gauss=zeros(sizex,sizey); %2D filter matrix
%gaussian filter
for i=-(sizex-1)/2:(sizex-1)/2
    for j=-(sizey-1)/2:(sizey-1)/2
        x0=(sizex+1)/2; %center
        y0=(sizey+1)/2; %center
        x=i+x0; %row
        y=j+y0; %col
        gauss(y,x)=exp(-(((x-x0).^2)/(2*sigma(1).^2)...
                     +((y-y0).^2)/(2*sigma(2).^2)...
                     ));   
    end
end
%normalize gaussian filter
gauss=gauss/sum(gauss(:));
end
%------------- END OF CODE --------------
