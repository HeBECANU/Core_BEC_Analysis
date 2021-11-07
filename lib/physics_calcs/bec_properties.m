function details=bec_properties(trap_freq,atom_number,mass,a_scat_len)
% calculate key properties of bose-einstein condensates
% input
%   omega-[scalar or 1x3] in Hz (not in rad/s)
%   atom number - scalar
% output
% 
% todo: Vectorize omega x atom_num

warning('use some caution with this function a full test suite is not yet developed')
%if the constants strucure already exists do not run
% set up the constants
global const
if isempty(const)
    hebec_constants
end

% input parsing
% omega
trap_freq=col_vec(trap_freq);
if numel(trap_freq)==1
    trap_freq=cat(1,trap_freq,trap_freq,trap_freq);
elseif numel(trap_freq)~=3
    error('omega must be 1, or 3 elements')
end

if sum(trap_freq<0)>0 && sum(imag(trap_freq)>0)>0
    error('omega must be positive real')
end

omega=2*pi*trap_freq;

%mass
if nargin<3 || isnan(mass) || isempty(mass)
    mass=const.mhe;
end

%scat_len
if nargin<4 || isnan(a_scat_len) || isempty(a_scat_len)
    a_scat_len=7.512000000000000e-09;
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

%pethick 5.54
u_eff_interaction=4*pi*(const.hb^2)*a_scat_len/mass;

n_max_peak_density=mu_chem_pot/u_eff_interaction;

%find the average density inside the elipsoid
tf_volume=(4/3)*pi*prod(r_tf_radi);
n_mean_density=atom_number/tf_volume;

tan_constant = (n_max_peak_density)*((64*pi^2)/7)*a_scat_len^2;

% goodness of TF p.171,2nd ed. pethick_smith
a_bar=sqrt(const.hb/(const.mhe*omega_bar)); % a bar is a the characteristic length
tf_goodness=atom_number*a_scat_len/a_bar;
% TF is good approx for bulk when this is much greater than 1

% the maximum velocity an atom placed at the top of the TF mean field potenial would reach in rolling off the potential
% this sets the bonds of the central disk of the PAL velocity distribution
% 1/2 m v^2 =U_tf
% v_max=sqrt(2*U_ft/m)
details.pal_vmax=sqrt(2*mu_chem_pot/mass);

% not sure which one is gravity axis so do all
grav_sag=-const.g0./(omega.^2);
details.tc.non_interacting=tc_non_interacting;
details.tc.finite_number=tc_finite_number;
details.tc.finite_interacting=tc_finite_interacting;
details.omega_mean=omega_mean;
details.omega_bar=omega_bar;
details.tf_radi=r_tf_radi;
details.tf_radi_bar=r_bar;
details.tf_radi_mean=r_mean;
details.tf_volume=tf_volume;
details.grav_sag=grav_sag;
details.a_bar=a_bar;
details.tf_goodness=tf_goodness;
details.density_peak=n_max_peak_density;
details.density_mean=n_mean_density;
details.u_eff_interaction=u_eff_interaction;
details.tan_contact= tan_constant;
details.inputs=[];
details.inputs.omega=omega;
details.inputs.atom_number=atom_number;
details.inputs.mass=mass;
details.inputs.a_scat_len=a_scat_len;



end

    

    
    
