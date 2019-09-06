function details=bec_properties(omega,atom_number,mass,a_scat_len)
% calculate key properties of bose-einstein condensates
% input
%   omega-[scalar or 1x3] in hz
%   atom number - scalar
% output
% 

warning('do not trust this function, it has not been tested well enough')
%if the constants strucure already exists do not run
% set up the constants
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

%scat_len
if nargin<4 || isnan(a_scat_len) || isempty(a_scat_len)
    a_scat_len=const.ahe_scat;
end


% begin calculation


details=[];

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
                      (1.33*a_scat_len*atom_number^(1/6)/(abar));

%fprintf('tc finite non interacting %g \n',tc_finite_interacting)                  


% other things we can calculate
% chemical potential
%pethick eq6.35
mu_chem_pot=((15^(2/5))/2)...
            *((atom_number*a_scat_len/abar)^(2/5))*...
            const.hb*omega_bar;
        
details.mu_chem_pot=mu_chem_pot;  

%pethick eq6.33     
r_tf_radi=sqrt(2*mu_chem_pot./(mass.*omega.^2));

r_bar=prod(r_tf_radi)^(1/3);
r_mean=mean(r_tf_radi);

%pethick eq6.33 
u_eff_interaction=4*pi*const.hb^2*a_scat_len/mass;

n_max_peak_density=mu_chem_pot/u_eff_interaction;

%find the average density inside the elipsoid
tf_volume=(4/3)*pi*prod(r_tf_radi);
n_mean_density=atom_number/tf_volume;


details.tc.non_interacting=tc_non_interacting;
details.tc.finite_number=tc_finite_number;
details.tc.finite_interacting=tc_finite_interacting;
details.omega_mean=omega_mean;
details.omega_bar=omega_bar;
details.tf_radi=r_tf_radi;
details.tf_radi_bar=r_bar;
details.tf_radi_mean=r_mean;
details.tf_volume=tf_volume;
details.density_peak=n_max_peak_density;
details.density_mean=n_mean_density;
details.u_eff_interaction=u_eff_interaction;



end

    

    
    
