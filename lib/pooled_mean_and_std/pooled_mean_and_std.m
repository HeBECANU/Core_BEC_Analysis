function [mean_out,std_out]=pooled_mean_and_std(mean_obs,std_obs,num_obs)
% pooled_mean_and_std - compute the standard deviation in a combned data set 
% from the standard deviations measured from data subsets
% uses expressions from 
% https://en.wikipedia.org/wiki/Pooled_variance
% https://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation
% this is an exact equivelence not an approximation

% Syntax:         pooled_mean_and_std(mean_obs,std_obs,num_obs)
% Inputs:
%    mean_obs   -  vector of the mean of each subset observation
%    std_obs    -  vector of the standard deviation in the subset observation
%    num_obs    -  vector of the population in each subset
%
% Outputs:
%    mean_out   - col vector of the mean in the combined data set
%    std_out    - col vector of the standard deviation in the combined data set
%
% Example: 
%     for example see test_pooled_mean_and_std
%
% Other m-files required: vol_vec
% Also See: test_pooled_mean_and_std
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%  - none so far
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2019-11-19

%------------- BEGIN CODE --------------

mean_obs=col_vec(mean_obs);
std_obs=col_vec(std_obs);
num_obs=col_vec(num_obs);


total_num_obs=sum(num_obs);

mean_out=sum(mean_obs.*num_obs)/total_num_obs;

% std_out=sqrt( (1/(total_num_obs-1))*...
%                 (   sum((num_obs-1).*(std_obs.^2) +  num_obs.*(mean_obs.^2) )  -...
%                     (mean_out.^2)*sum(num_obs) ...
%                 )...
%              )
% an equivelent way that is a little easier to understand
std_out=sqrt( (1/(total_num_obs-1))*...
                (   sum((num_obs-1).*(std_obs.^2) +  num_obs.*((mean_obs-mean_out).^2) ) ...
                )...
             );


end