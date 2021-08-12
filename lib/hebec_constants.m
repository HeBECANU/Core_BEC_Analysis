%constants
%TODO
% - add subsections eg SI,imperial, customary
% - include the units
% - make this a class so the values cant be set
% - add references for all values
% - include transitions

function const = hebec_constants()
% global const
%fundamental
const.c = 299792458; %speed of light (m/s)
const.h = 6.626070040*10^-34; %Planck constant (J s)
const.hb = 1.054*10^-34; %reduced Planck constant (J S)
const.kb = 1.3806488*10^-23; %Boltzmann constant (*m2 kg s-2 K-1*)
const.mu0 = 1.2566370614*10^-6;  % vacuum permeability [Tm/A]
const.epsilon0 = 8.858*10^-12;%(*electric permittivity of free space*)
%elemental
const.mub =9.274009994*10^-24; %Bohr magneton*(J/T)
const.electron = 1.60217657*10^-19;%(*charge of electron*)
const.me = 9.10938291*10^-31;%(*mass of electron*)
const.mp = 1.67262158*10^-27;%(*mass of proton*)
const.grav=6.67430*10^-11;  %Newtonian constant of gravitation %https://physics.nist.gov/cgi-bin/cuu/Value?bg

%Helium
% const.ahe_scat=15*10^-9;
const.ahe_scat=7.512000000000000e-09; %m^2 Moal et al https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.023203
const.b_freq=2*const.mub/const.h /(1e4); %for the metastable state - hz/G
% calculated for the metastable state
const.mhe = 1.66*10^-27*4.002;%(*helium mass*)
const.interaction_energy = 4*pi*const.hb^2*const.ahe_scat/const.mhe; % interaction strength []

%environmental
const.g0 =9.7960207-78e-8;
% (-35 17.0', 149 6.8) alt=560.78m
% https://d28rz98at9flks.cloudfront.net/15167/Rep_261.pdf 
%with correction of 78e-8 from http://www.publish.csiro.au/ex/ASEG2007ab149 applied
%9.796 cbr approx grav at (-35.283103 , 149.112634) from http://ddfe.curtin.edu.au/gravitymodels/GGMplus/data/ga/S40E145.ga.png

%customary
const.a0 = 0.529*10^-10;%(*bohr radius*)

%experiment
const.fall_distance = 0.8587; 

const.mu = 9.27e-28; %J/G
const.h = 6.62607015e-34;
const.hbar = const.h/(2*pi);
const.f_mu = const.mu/const.h;
const.w_mu = const.mu/const.hbar;
const.c = 299792458;
const.q = 1.602e-19;
% Notation & lookup
const.terms = {'S','P','D','F','G','H','I','K'};
%% Reference values
const.f_table.g_2_3P_2.e_5_3S_1 = 727.3032446e12;
% Misc transitions - what do the stars mean?
const.f_table.g_2_3P_2.e_5_3P_0 = 1e9*const.c/404.628937550957;
const.f_table.g_2_3P_2.e_5_3P_1 = 1e9*const.c/404.629844755577;
const.f_table.g_2_3P_2.e_5_3P_2 = 1e9*const.c/404.629918705477; 
% Historically controversial transitions
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620968e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622889e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651246e12; 
% Singlet-triplet transitions
const.f_table.g_2_3P_2.e_5_1S_0 = 1e9*const.c/406.8886971706;
const.f_table.g_2_3P_2.e_5_1P_1 = 1e9*const.c/402.322271224483;
const.f_table.g_2_3P_2.e_5_1D_2 = 744.43034335e12; % 402.7nm

%Fitted valuse for the 5^3D's
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620836e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622758e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651114e12;

end
