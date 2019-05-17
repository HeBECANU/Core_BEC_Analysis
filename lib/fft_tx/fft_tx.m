function[out]=fft_tx(t,x,varargin)
% fft_tx - general purpose fft
% built to hande most use cases: uneven data, windowing & padding to increase the freq sampling. 
% all built to keep the returned amplitudes correct, this will mean that the integerated power 
% is wrong.
% https://dsp.stackexchange.com/questions/7788/setup-frequency-array-properly
% Uses the spread in the sampling times to dynamicaly change the resampling
% ratio

% Syntax:   freq_amp=fft_tx(times,val)
%           freq_amp=fft_tx(times,val,'padding',10,'window','gauss','win_param',{5});
%
% Inputs:
%   t          - vector, sample times (will deal with unmatched vector dimension)
%   x          - vector, the signal   (will deal with unmatched vector dimension)
%   Optional Name Value Pairs:
%       'padding'        - float, value>=1, factor to padd the data by to increase the frequeny sampling (but not resolvability)
%       'window'         - string, windowing function to apply to the data before processing.
%                       ['none','hamming','gauss','blackman','hanning']
%       win_param       - cell array, extra inputs after the length argument to the windowing function
%                           - gauss ,win_param={3}, reciprocal of the standard deviation
%                           - all others (besides none) 'symmetric'(recomended) or 'periodic'

% Outputs:
%    freq_amp - 2*L matrix
%               - freq_amp(1,:) are the frequency bins
%               - freq_amp(2,:) are the (complex) amplitudes
%
% Example: 
%         times=linspace(0,1e3,1e6);
%         testfun=@(t) 100*sin(2*pi*t*100+pi)+1*sin(2*pi*t*133+pi);
%         val=testfun(times);
%         out=fft_tx(times,val,'padding',10,'window','gauss','win_param',{5});
%         figure(5)
%         plot(out(1,:),abs(out(2,:)))

% Other m-files required: none
% Also See: test_fft_tx
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%    - in out as struct
%    - add power and amplitude normalization with tests
%    - better default window params
%    - skipable sort for inputs that are already sorted
%    - skipable regular sampling check
%    - chosable resampling stratagey
%    - add verbosity option
%    - use x and t as col vectors in the code
%    - add option to padd to power of 2 for more speed if its close
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-11-21

%------------- BEGIN CODE --------------

%% turn x,t inputs into row vectors
%adaptively deal with the data if its in row or col fromat
%this method of checking sizes is faster than doing x(:)
% - test with x,t 2e6 elements
%   - input wrong way (col vec)
%     - case method ~14ms
%     - x(:) method ~17ms
%   - input right way (row vec)
%     - case method ~100us
%     - x(:) method ~18ms
% TODO change code to use col vectors

if size(size(t),2)==2
    if size(t,1)~=1 && size(t,2)==1
        t=t';
    elseif size(t,1)<=1 && size(t,2)<=1
        error('thats not a vector in t')
    end
else
    error('you have tried to input the wrong shape in t')
end
if size(size(x),2)==2
    if size(x,1)~=1 && size(x,2)==1
        x=x';
    elseif size(x,1)==1 && size(x,2)==1
        error('thats not a vector in x')
    end
else
    error('you have tried to input the wrong shape in x')
end  


%% parse the optional inputs
win_ok_fun=@(x) sum(contains({'none','hamming','gauss','blackman','hanning','kaiser','chebyshev','bohmanwin','blackmanharris'},x))==1;
is_c_logical=@(in) isequal(in,true) || isequal(in,false); %can x be cast as a logical
p = inputParser;
addParameter(p,'window','none',win_ok_fun);
addParameter(p,'win_param',{},@(x) true);
addParameter(p,'padding','',@(x) isfloat(x) && x>=1);
addParameter(p,'issorted',0,is_c_logical);
addParameter(p,'verbose',1,@(x) round(x)==x && x>=0);
parse(p,varargin{:});
pad=p.Results.padding;
window_fun=p.Results.window;
window_param=p.Results.win_param;
xsorted=p.Results.issorted;
verbose=p.Results.verbose;

%% sort the data if needed
if ~xsorted
    [t,sort_idx]=sort(t);
    x=x(sort_idx);
end

%% resample the data if needed
% TODO: some user options for resampling stratagies

% find if the range of subsequent sample time differences exceeds machine error
% linspace generated tvecs have a range of 2*eps(tend)
time_diffs=diff(t); 
if range(time_diffs)>10*eps(t(end)) 
    if verbose>=1
        fprintf('%s:resampling inputs \n',mfilename);
    end
    %TODO: this algorithm needs improving
    %  https://escholarship.org/uc/item/4rb242mv
    %  https://www.math.ucdavis.edu/~strohmer/research/sampling/irsampl.html
    %  https://www.eurasip.org/Proceedings/Eusipco/Eusipco2011/papers/1569428115.pdf
    dyn_res_factor=10+10*abs(std(time_diffs)/mean(time_diffs));
    num_resamp=2^nextpow2(numel(t)*dyn_res_factor); %for faster fft wa can round to powers of 2
    t_resample=linspace(min(t),max(t),num_resamp);
    x=interp1(t,x,t_resample,'spline');
    t=t_resample;
end


dt = (t(2)-t(1));             % Sampling period
%disp(num2str(dt))
fs=1/dt;
len_before_pad = numel(x);             % Length of signal
%apply windowing function
switch window_fun
    case 'hamming'
        win=hamming(len_before_pad,window_param{:})';
    case 'gauss'
        win=gausswin(len_before_pad,window_param{:})'; 
    case 'blackman'
        win=blackman(len_before_pad,window_param{:})'; 
    case 'hanning'
        win=hann(len_before_pad,window_param{:})'; 
    case 'kaiser'
        win=kaiser(len_before_pad,window_param{:})';
    case 'chebyshev'
        win=chebwin(len_before_pad,window_param{:})';
    case 'bohmanwin'
        win=bohmanwin(len_before_pad)';
    case 'blackmanharris'
        win=blackmanharris(len_before_pad)';
    case 'none'
        %do nothing
end
if ~isequal(window_fun,'none')
    win=win.*(len_before_pad/sum(win));
    x=x.*win;
end

if pad<1
    error('pad cant be less than one')
elseif pad~=1
    x=[x,zeros(1,round(len_before_pad*(pad-1)))];
end

len = numel(x);  % Length of signal
%disp(num2str(L))
%t = (0:L-1)*T;        % Time vector
y = fft(x);
amp = y/len_before_pad;
if  mod(len,2) %odd case
    niq=ceil(len/2);
    amp = amp(1:niq);
    amp(2:end) = 2*amp(2:end); %multiply all non DC values by 2 because of half cut
    f = fs/len*((0:niq-1));
else
    niq=len/2+1;
    amp = amp(1:niq);
    amp(2:end) = 2*amp(2:end); %multiply all non DC values by 2 because of half cut
    f = fs*(0:(len/2))/len;  
end


out=[f;amp];


end


%% alternate way of turning inputs into row vec
% if sum(size(x)==1)==0
%  error('you have tried to input the wrong shape in x')
% end
% if sum(size(t)==1)==0
%  error('you have tried to input the wrong shape in t')
% end
% x=x(:)';
% t=t(:)';
