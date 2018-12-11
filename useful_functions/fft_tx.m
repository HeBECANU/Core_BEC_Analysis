function[out]=fft_tx(t,x)
%a fft function that can handle unevenly spaced data
%uses the spread in the sampling times to dynamicaly change the resampling
%ratio
%https://dsp.stackexchange.com/questions/7788/setup-frequency-array-properly


%first sort the data
[t,i]=sort(t);
x=x(i);
%find the number of sampling times present
sample_times=uniquetol(t(1:end-1)-t(2:end),1e-10); %avoid machine error
if size(sample_times,1)>1 %if not uniform
    fprintf('resampling\n');
    dyn_res_factor=10+500*abs(std(sample_times)/mean(sample_times));
    t_resample=linspace(min(t),max(t),size(t,1)*dyn_res_factor);
    x=interp1(t,x,t_resample,'spline')';
    t=t_resample;
end
dt = (t(2)-t(1));             % Sampling period
%disp(num2str(dt))
fs=1/dt;
len = size(x,1);             % Length of signal
%disp(num2str(L))
%t = (0:L-1)*T;        % Time vector
y = fft(x);
p2 = abs(y/len);

if  mod(len,2) %odd case
    niq=ceil(len/2);
    p1 = p2(1:niq);
    p1(2:end) = 2*p1(2:end);
    f = fs/len*((0:niq-1));
else
    niq=len/2+1;
    p1 = p2(1:niq);
    p1(2:end) = 2*p1(2:end); %multiply all non DC values by 2 because of half cut
    f = fs*(0:(len/2))/len;  
end

f=transpose(f);
out=[f p1];


end
