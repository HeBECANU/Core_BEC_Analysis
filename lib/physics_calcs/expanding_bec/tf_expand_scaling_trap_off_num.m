function [lambda_end,soln_sequence,lambda_fun]=tf_expand_scaling_trap_off_num(omega_tzero,tmax)
% calculates the rescaling paramters that are given by the ODE
%  \derivn[2]{\lambda_j}{t}= \frac{\omega_j^2(0)}{\lambda_j \lambda_1 \lambda_2 \lambda_3} - \omega_j^2(t) \lambda_j
% in the case that  \omega_j(t>0)=0
%  \derivn[2]{\lambda_j}{t}= \frac{\omega_j^2(0)}{\lambda_j \lambda_1 \lambda_2 \lambda_3}
% we want the solution of lambda from t=0 to t=tmax
omega_tzero=col_vec(omega_tzero);
omega_max=max(omega_tzero);




%ode45
%ode113
inital_lambda_dlambda=[1,1,1,0,0,0]';
twindow=[0,tmax];

% ode_opts = odeset('RelTol',1e-9,'AbsTol',1e-10,'Stats','on','InitialStep',1e-3/omega_max,'MaxStep',10/omega_max,'Vectorized',0); %
% [t_soln,lambda_dlambda_soln] = ode45(@(t,X) rescaling_deriv(t,X,omega_tzero),twindow,inital_lambda_dlambda,ode_opts);


ode_opts = odeset('RelTol',1e-9,'AbsTol',1e-12,'Stats','off','InitialStep',1e-3/omega_max,'MaxStep',10/omega_max,'Vectorized',0); %
lambda_dlambda_soln = ode113(@(t,X) rescaling_deriv(t,X,omega_tzero),twindow,inital_lambda_dlambda,ode_opts);

%txvdata = [t,y(:,1),y(:,2)];
if lambda_dlambda_soln.x(end)~=tmax
    error('did not evaluate at this point')
end
lambda_end=lambda_dlambda_soln.y(end,1:3);

soln_sequence=[];
soln_sequence.time=lambda_dlambda_soln.x;
soln_sequence.lambda=lambda_dlambda_soln.y(:,1:3);
soln_sequence.d_lambda=lambda_dlambda_soln.y(end,4:6);

% for handy use later we pass out a function which can be used to evaluate the function at any time
if nargout>2
    lambda_fun=@(t) evaluating_fun(lambda_dlambda_soln,t);
end

end



function dL1_dL2 = rescaling_deriv(t,L1_L2,omega_tzero)
    % this function is used in the ode solver to get a solution for
    %  \derivn[2]{\lambda_j}{t}= \frac{\omega_j^2(0)}{\lambda_j \lambda_1 \lambda_2 \lambda_3}
    % to solve this ode we use a change of variable
    % L1=lambda          (where L1 is a 1x3 vector)
    % L2= d lambda/ d t
    % therfore
    % d L1 /d t = L2
    % d L2/ d t = \frac{\omega_j^2(0)}{L1(j) L1(1) L1(2) L1(3)}
    L2=L1_L2(4:end);
    L1=L1_L2(1:3);
    
    dL1=L2;
    dL2= ((omega_tzero.^2)./L1) *(1/(prod(L1)));
    
    dL1_dL2=cat(1,dL1,dL2);

end %function


function lambda=evaluating_fun(lambda_dlambda_soln,t)
    t=col_vec(t);
    % used to extrapolate outside the simulation interval
    tmin=lambda_dlambda_soln.x(1);
    tmax=lambda_dlambda_soln.x(end);
    mask_before_sim= t<tmin;
    mask_after_sim=t>tmax;
    mask_in_sim=~mask_before_sim & ~mask_after_sim;
    lambda=nan(numel(t),3);
    % evaluate the times in the sim
    if sum(mask_in_sim)>0
        lambda_dlambda_tmp=deval(lambda_dlambda_soln,t(mask_in_sim));
        lambda(mask_in_sim,:)=transpose(lambda_dlambda_tmp(1:3,:));
    end
    if sum(mask_before_sim)>0
         lambda(mask_before_sim,:)=repmat(transpose(lambda_dlambda_soln.y(1:3,1)),sum(mask_before_sim),1);
    end
    if sum(mask_after_sim)>0
        % if adter sim we use a very simple linear interp
        % L=L_end+(t-tend)*dL/dt|_t=tend
        dt=t(mask_after_sim)-tmax;
        lambda(mask_before_sim,:)=repmat(transpose(lambda_dlambda_soln.y(1:3,end)),sum(mask_before_sim),1)+...
            repmat(dt,1,3).*repmat(transpose(lambda_dlambda_soln.y(4:6,end)),sum(mask_before_sim),1);
    end
    

end