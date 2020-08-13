function [alpha_au,alpha_si]=aprox_he_polz(optical_freq)
% find the polarizability of He* using published osc strengths
% built for the TO project meant more for finding the curvature arround the TO that the intercept itself
% example
% aprox_he_polz(linspace(1,600,1e2)*1e12)
% TODO:
% - use drake tables for states above n=3
% - add inital and final state information to the metastable_helium_transition_data structure
%       having the qunatum numbers for each state would be handy
% example
%   optical_freq=linspace(1,600,1e2)*1e12;
% Bryce Henson 2020-07-09
%%


optical_freq=col_vec(optical_freq);
global const

atom_u=[];
% atom_u.energy=(const.hb)^2/(const.me*(const.a0^2)); 
% atom_u.time=const.hb/atom_u.energy;
atom_u.energy=4.3597447222071e-18 ;%J
atom_u.time=2.4188843265857e-17 ;%s
atom_u.polz=1.64877727436e-41;% C2⋅m2⋅J−1 
atom_u.second_hyp_polz=6.2353799905e-65;% C^4⋅m^4⋅J^-3
atom_u.a0=5.29177210903e-11 ;% m

% use the theory calculated value of the polz for states above n=3
% found by comparing simple model with the 315.54238 value in https://journals.aps.org/pra/pdf/10.1103/PhysRevA.93.052516
% this means that the polarizability will not be exactly correct
%n_gt_3_polz=4.844735; %a_0^2 cgs vol polz units
n_gt_3_polz=7.065; %a_0^2 cgs vol polz units
transitions=metastable_helium_transition_data;


% first we preprocess the structure to get an array of osc strength and transition freq
%dipole_mask=cellfun(@(x) strcmp(x.type,'E1'),transitions);
dipole_mask=true(numel(transitions),1);

osc_st_freq_arr=cat(2,  cellfun(@(x) x.osc_strength.val, transitions(dipole_mask)),...
                    cellfun(@(x) x.frequency_hz.val, transitions(dipole_mask)) );
                           

freq_sq_diff_matrix=bsxfun(@minus, (osc_st_freq_arr(:,2)').^2, optical_freq.^2);

sq_energy_diff_atomic_u=freq_sq_diff_matrix*(const.h/atom_u.energy)^2;


polz_components=repmat(osc_st_freq_arr(:,1)',[size(optical_freq,1),1])./sq_energy_diff_atomic_u;
polz=sum(polz_components,2)+n_gt_3_polz;
alpha_au=    polz  ; % static_polz
alpha_si=alpha_au*atom_u.polz; 

                

