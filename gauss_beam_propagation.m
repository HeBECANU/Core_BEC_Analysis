wavelen = 1083.3e-9;

% starting z position (waist location) for x profile
zy=0.00000001;%0.141;%

% initial x waist and radius of curvature of the beam at focus, in m 
w0y=109e-6;%178e-6;%
zry = Rayleigh_len(w0y,wavelen);
R0y=zy*(1+(zry/zy)^2);

%desired beam siez
w2y = 1.42e-3;%1000e-6;%

% Raylenx=Rayleigh_len(w0x,lambda_dip)
% divx=divergence(lambda_dip,w0x)


f1 = 0.05;
f2 = 0.1;
z0 = (pi*w2y*w0y/wavelen*f1/f2+f1);%8.133743243150532e-01;%300e-3+5.95e-1;
dist1 = f1+f2+f1.^2/(z0-f1);

% z0 = z0-0.141;

fig_num = 101;

figure(fig_num)
clf

[R,w,z]=free_prop_gauss(w0y,R0y,zy,z0,wavelen,fig_num);
[R,w,z]=lens_gauss(w,R,z,f1,wavelen,fig_num);
[R,w,z]=free_prop_gauss(w,R,z,dist1,wavelen,fig_num);
[R,w,z]=lens_gauss(w,R,z,f2,wavelen,fig_num);
[R,w,z]=free_prop_gauss(w,R,z,1.245-z,wavelen,fig_num);

% ryleigh length
function Rz = Rayleigh_len(w0,wavelen)
Rz = pi*w0^2/wavelen;
end

% complex gaussian beam parameter
function out = R_fq(q)
    out= 1./(real(1./q));
end

% w complex gaussian beam parameter
function out = w_fq(q,wavelen)
    out = sqrt(wavelen./(pi.*imag(-1./q)));
end

% complex gaussian beam parameter
function q = q_param(R,w,wavelen)
    q= 1/(1/R-1j*wavelen/(pi*w^2));
end

function [R_new,w_new,z] = free_prop_gauss(w,R,z0,L,wavelen,fig_num)
    q0=q_param(R,w,wavelen);
    
    L_step=L/1000.0;
    L_arr=linspace(0.0,L+L_step,round((L+L_step)/L_step));
    z_arr=L_arr+z0;
    
    z=z0+L;
    
    q_arr=q0+L_arr;    
    R_arr = R_fq(q_arr);
    w_arr=w_fq(q_arr,wavelen);
    
    figure(fig_num)
    hold on
    plot(z_arr,w_arr,'r-')
    hold off
    
    R_new = R_arr(end);
    w_new = w_arr(end);
%     out= [R_arr(end),w_arr(end),z];
end

function [R,w,z] = lens_gauss(w,R,z,f,wavelen,fig_num)
    q0=q_param(R,w,wavelen);
    
    q=1/(1/q0-1/f);
    R= R_fq(q);
    w=w_fq(q,wavelen);
    
    figure(fig_num)
    hold on
%     plt.axvline(x=z,ymin=0,ymax=1, linewidth=5, color='b',zorder=1)
    plot([z,z],[0,w],'b-')
    hold off
%     out = ;
    end