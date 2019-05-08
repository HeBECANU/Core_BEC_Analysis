%there are lots of windowing functions that the user can call
% to improve usability we want the defaults to be fairly similar
% this script tries to find defaults for asjustable windows so that they are similar to gausswin(len,5)
% it would be great to come up with a window converter that converts some gaussian to any number of window function
% parameters
%% demonstration of the problem
len=1e3;
plot(hamming(len))
hold on
plot(gausswin(len)) %adjustable
plot(blackman(len))
plot(hann(len))
plot(kaiser(len)) %adjustable
plot(chebwin(len)) %adjustable
%plot(bohmanwin(len)) 
plot(blackmanharris(len))
hold off 
legend('hamming','gausswin','blackman','hann','kaiser','chebwin','blackmanharris')

%%


wvtool(hamming(len),...
gausswin(len,5),... %adjustable
blackman(len),...
hann(len),...
kaiser(len,4),... %adjustable
chebwin(len,200),...
blackmanharris(len)) %adjustable

%% manual adj

figure(1)
clf
len=1e3;
plot(hamming(len))
hold on
plot(gausswin(len,3)) %adjustable
plot(blackman(len))
plot(hann(len))
plot(kaiser(len,7)) %adjustable
plot(chebwin(len,100)) %adjustable
hold off
legend('hamming','gausswin','blackman','hann','kaiser','chebwin')


%%
gwidth=4;
taylor_numlobes=4;
param_cheb=match_window_funs(@(n) gausswin(n,gwidth), @(n,p) chebwin(n,p),1e3,100,2,0);
param_kaiser=match_window_funs(@(n) gausswin(n,gwidth), @(n,p) kaiser(n,p),1e3,3,2,0);

figure(1);
clf
len=1e3;
plot(gausswin(len,gwidth)) %adjustable
hold on
plot(kaiser(len,param_kaiser)) %adjustable
plot(chebwin(len,param_cheb)) %adjustable
plot(hanning(len))
plot(hamming(len))
hold off
legend('gausswin','kaiser','chebwin','hanning','hamming')

wvtool(gausswin(len,gwidth),... %adjustable
kaiser(len,param_kaiser),... %adjustable
chebwin(len,param_cheb))%,... %adjustable
%hanning(len),...
%hamming(len))


%% convert the non adjustable windows into the closest gaussian window

match_window_funs(@(n) hamming(len), @(n,p) gausswin(n,p),1e3,2,2,2)

%%

match_window_funs(@(n) blackman(n), @(n,p) gausswin(n,p),1e3,3,2,2)

%%
match_window_funs(@(n) hanning(n), @(n,p) gausswin(n,p),1e3,3,2,2)

%%

match_window_funs(@(n) blackmanharris(n), @(n,p) gausswin(n,p),1e3,3,2,2)




