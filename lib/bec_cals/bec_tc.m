function [tc_out,details]=bec_tc(omega,atom_number,mass)
% input
%   omega-[scalar or 1x3] in hz
warning('this function will be merged into bec_properties, please move your code use')

%if the constants strucure already exists do not run
global const
if isempty(const)
    hebec_constants
end

% input parsing
% omega
omega=col_vec(omega);
if numel(omega)==1
    omega=cat(1,omega,omega,omega);
elseif numel(omega)~=3
    error('omega must be 1, or 3 elements')
end

if sum(omega<0)>0 && sum(imag(omega)>0)>0
    error('omega must be positive real')
end

omega=2*pi*omega;

%mass
if nargin<3 || isnan(mass) || isempty(mass)
    mass=const.mhe;
end



omega_bar=prod(omega)^(1/3);
omega_mean=mean(omega);


tc_non_interacting=const.hb *omega_bar*(atom_number^(1/3))/((zeta(3)^(1/3))*const.kb);

%fprintf('tc non interaging %g \n',tc_non_interacting)
%pethick eq2 .91
tc_finite_number=tc_non_interacting*(1-0.73*omega_mean/((atom_number^(1/3))*omega_bar));

%fprintf('tc finite number %g \n',tc_finite_number)
abar=sqrt(const.hb/(mass*omega_bar));
%pethick eq11.14
tc_finite_interacting=tc_finite_number-tc_non_interacting*...
                      (1.33*const.ahe_scat*atom_number^(1/6)/(abar));

%fprintf('tc finite non interacting %g \n',tc_finite_interacting)                  
                  
tc_out=tc_finite_interacting;

details.tc_non_interacting=tc_non_interacting;
details.tc_finite_number=tc_finite_number;
details.tc_finite_interacting=tc_finite_interacting;


end

    

    
    
