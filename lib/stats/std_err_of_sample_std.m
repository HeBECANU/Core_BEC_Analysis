function out=std_err_of_sample_std(n,stdin)
%standard error of the sample standard deviation
%https://stats.stackexchange.com/questions/631/standard-deviation-of-standard-deviation
%https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation
% there are also some notes here http://davidmlane.com/hyperstat/A19196.html
% https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
%gamma((n-1)/2)
%https://math.stackexchange.com/questions/1259383/calculating-uncertainty-in-standard-deviation

%warning('needs testing, will be biased for small N')
out=stdin/sqrt((2*n-2));


end
