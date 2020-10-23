centre_val = 0.1;
samp = randn(1e5,1)+centre_val;
edges = linspace(-3,3,201);
C = histcounts(samp,edges);
T = midpoints(edges);
V = diff(edges);
Y = C./V;
% make saturated model...

thr = 0.1;
% Ycap(rescale(Y)>thr) = 
m = rescale(Y)>thr;
Ydash = Y;
Ydash(m) = thr*max(Y);
idx = mean(T(m));
val = mean(T.*Ydash)/sum(Ydash);

A = .4e-6;
tau = 2.;
quash_fun = @(t) [zeros(size(t)),A*exp(-t/tau)];
tdash = 0:mean(diff(T)):10;
qe_quash = 1-conv(Y,quash_fun(tdash),'same');
qe_quash(qe_quash<0) = 0;

Y_sat = qe_quash.*Y;
[v_sat,i_sat] = sat_pulse_centre(T,Y_sat,0.1);

m_sat = rescale(Y_sat)>thr;
Ydash_sat = Y_sat;
Ydash_sat(m_sat) = thr*max(Y_sat);
idx_sat = mean(T(m_sat));
val_sat = mean(T.*Ydash_sat)/sum(Ydash_sat);

cli_header('Saturated peak centering:');
cli_header(1,'Support COM error: %.3e',idx_sat-centre_val);
cli_header(1,'Threshold COM error: %.3e',val_sat-centre_val);

stfig('Saturated centering tests');
clf
subplot(2,2,1)
hold on
plot(T,Y)
plot(T,Ydash)
plot(T(m),thr*max(Y)*ones(size(T(m))),'k')
plot(idx,0,'rx')
plot(val,0,'ko')
title('Raw profile')

ylabel('Density')
xlabel('T')
subplot(2,2,2)
plot(quash_fun(tdash))
title('Accumulation function')

subplot(2,2,3)
plot(T,qe_quash)
title('QE loss envelope')

subplot(2,2,4)
hold on
plot(T,Y_sat)
plot(T,Y,'k:')
plot(T,Ydash_sat)
plot(T(m_sat),thr*max(Y_sat)*ones(size(T(m_sat))),'k')
plot(T(i_sat),0,'rx')
plot(v_sat,0,'ko')
xlim([-5,5])
title('Saturated peak')

% first: threshold, find COM

    