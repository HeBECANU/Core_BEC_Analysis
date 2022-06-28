f1 = 200;%focal length of first lens
f2 = 125;%focal length of second lens

do1 = 205+180;%397;%375;%mm
d_lens = 230+390;%455+50;%490:5:550;%450+105-30;%570;

di1 = 1./(1./f1-1./do1);

do2 = d_lens - di1;
di2 = 1./(1./f2-1./do2)

m1 = -di1./do1;
m2 = -di2./do2;
m = m1*m2