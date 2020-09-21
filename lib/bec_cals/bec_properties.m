function details=bec_properties(omega,atom_number,varargin)
% calculate key properties of bose-einstein condensates
% input
%   omega-[scalar or 1x3] in hz
%   atom number - scalar
% optional keyword arguments
%       mass - particle mass (default const.mhe = 6.6433e-27 kg
%       a_scat_len - s-wave scattering length, default 7.512e-9 m
%       temperature - Default 0 K
% output
% 
% todo: Vectorize omega x atom_num

warning('do not trust this function, it has not been tested well enough')
warning('bec_properties now accepts inputs (trap_freqs, total_number,varargin')

%if the constants strucure already exists do not run
% set up the constants
global const
if isempty(const)
    const = hebec_constants;
end

p = inputParser;
%     defaultFontSize = 12; 
%     defaultFont = 'times'; 
addParameter(p,'temperature',0);
addParameter(p,'mass',const.mhe);
addParameter(p,'a_scat_len',const.ahe_scat);
parse(p,varargin{:});
temperature = p.Results.temperature;
mass = p.Results.mass;
a_scat_len = p.Results.a_scat_len;

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


% begin calculation
details=[];

omega_bar=prod(omega)^(1/3);
omega_mean=mean(omega);
atom_number = col_vec(atom_number);
zeta_3 = 1.202056903;
tc_non_interacting=const.hb *omega_bar*(atom_number.^(1/3))/((zeta_3^(1/3))*const.kb);

%fprintf('tc non interaging %g \n',tc_non_interacting)
%pethick eq2 .91
tc_finite_number=tc_non_interacting.*(1-0.73*omega_mean./((atom_number.^(1/3))*omega_bar));

%fprintf('tc finite number %g \n',tc_finite_number)
abar=sqrt(const.hb/(mass*omega_bar));
%pethick eq11.14
tc_finite_interacting=tc_finite_number-tc_non_interacting.*...
                      (1.33*a_scat_len*atom_number.^(1/6)/(abar));


%fprintf('tc finite non interacting %g \n',tc_finite_interacting)                  

condensed_fraction = 1 - (temperature./tc_finite_interacting).^3;
condensed_number = atom_number.*condensed_fraction; 

% other things we can calculate
% chemical potential
%pethick eq6.35
mu_chem_pot=((15^(2/5))/2)...
            *((condensed_number *a_scat_len/abar).^(2/5))*...
            const.hb*omega_bar;
        
details.mu_chem_pot=mu_chem_pot;  

%pethick eq6.33     
r_tf_radi=sqrt(2*mu_chem_pot./(mass.*omega'.^2));

r_bar=prod(r_tf_radi).^(1/3);
r_mean=mean(r_tf_radi);

%pethick eq6.33 
u_eff_interaction=4*pi*const.hb^2*a_scat_len/mass;

n_peak_density=mu_chem_pot/u_eff_interaction; % 

%find the average density inside the elipsoid
tf_volume=(4/3)*pi*prod(r_tf_radi,2);
n_mean_density=condensed_number ./tf_volume; %Pethick % Smith 6.31

tan_contact= condensed_number .*(n_peak_density)*((64*pi^2)/7)*a_scat_len^2; %Chang 2016 & Ross 2020(ish)

LHY_correction.mean = 2.5*(128/(15*sqrt(pi)))*sqrt(n_mean_density*a_scat_len^3);
LHY_correction.peak = 2.5*(128/(15*sqrt(pi)))*sqrt(n_peak_density*a_scat_len^3);


speed_of_sound = sqrt(u_eff_interaction*n_peak_density/mass); %Pethick and Smith 8.92
oscillator_length = sqrt(const.hb./(mass*omega));
healing_length = 1/(8*pi*a_scat_len*n_peak_density); % Pethick % Smith 6.62

% landau_critical_velocity = 
% Baym & Pethick 2012 https://arxiv.org/pdf/1206.7066.pdf

details.tc.non_interacting=tc_non_interacting;
details.tc.finite_number=tc_finite_number;
details.tc.finite_interacting=tc_finite_interacting;
details.omega_mean=omega_mean;
details.omega_bar=omega_bar;
details.tf_radi=r_tf_radi;
details.tf_radi_bar=r_bar;
details.tf_radi_mean=r_mean;
details.tf_volume=tf_volume;
details.density_peak=n_peak_density;
details.density_mean=n_mean_density;
details.u_eff_interaction=u_eff_interaction;
details.tan_contact= tan_contact;
details.LHY_correction= LHY_correction;
% details.landau_limit = landau_limit;
details.healing_length = healing_length;
details.oscillator_length = oscillator_length;
details.speed_of_sound = speed_of_sound;
details.condensed_fraction = condensed_fraction;
details.lambda_db = deBroglie_wavelength(temperature);

end

    

    
    
