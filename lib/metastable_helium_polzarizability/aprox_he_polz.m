function [alpha_au,alpha_si,higher_terms]=aprox_he_polz(optical_freq,transitions)
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
global const

if nargin<2
    transitions=metastable_helium_transition_data;
end

%apply_polz_offset=true;
apply_polz_offset=false;

optical_freq=col_vec(optical_freq);

atom_u=atomic_units;

% use the theory calculated value of the polz for states above n=3
% found by comparing simple model with the 315.54238 value in https://journals.aps.org/pra/pdf/10.1103/PhysRevA.93.052516
% this means that the polarizability will not be exactly correct
%n_gt_3_polz=4.844735; %a_0^2 cgs vol polz units
n_gt_3_polz=7.065; %a_0^2 cgs vol polz units



% first we preprocess the structure to get an array of osc strength and transition freq
%dipole_mask=cellfun(@(x) strcmp(x.type,'E1'),transitions);
dipole_mask=true(numel(transitions),1);

osc_st_freq_arr=cat(2,  cellfun(@(x) x.osc_strength.val, transitions(dipole_mask)),...
                    cellfun(@(x) x.frequency_hz.val, transitions(dipole_mask)) );
                           

freq_sq_diff_matrix=bsxfun(@minus, (osc_st_freq_arr(:,2)').^2, optical_freq.^2);

sq_energy_diff_atomic_u=freq_sq_diff_matrix*(const.h/atom_u.energy)^2;


polz_components=repmat(osc_st_freq_arr(:,1)',[size(optical_freq,1),1])./sq_energy_diff_atomic_u;
polz=sum(polz_components,2);
if apply_polz_offset
    polz=polz+n_gt_3_polz;
end
alpha_au=    polz  ; % static_polz
alpha_si=alpha_au*atom_u.polz; 

%% vector and tensor terms
% to compute the vector and tensor we follow
% https://doi.org/10.1140/epjd/e2013-30729-x

% first we must convert the osc strengths into reduced matrix elements
osc_st_vec= cellfun(@(x) x.osc_strength.val, transitions(dipole_mask));
% need to use the quantum number of the 
freq_trans_vec= cellfun(@(x) x.frequency_hz.val, transitions(dipole_mask));
energy_diff_atomic_u=freq_sq_diff_matrix*(const.h/atom_u.energy);

qnum_ground_j=1;
qnum_excited_j_vec=cellfun(@(x) x.name(2+inlidx(regexp(x.name,'_\{\d{1}\}'),2)) , transitions(dipole_mask));
qnum_excited_j_vec=str2num(qnum_excited_j_vec);

% equation 13
red_mat_elm_sq = osc_st_vec ...
            .*(2.*qnum_ground_j+1)  ...
            .* (3*const.hb*(const.electron^2))./(2*const.me*(2*pi*freq_trans_vec));

% equation 11
polz_k_array=nan(numel(optical_freq),3);
iimax=numel(transitions);
for K=0:2
    prefactor= (-1)^(K+qnum_ground_j+1)*sqrt(2*K+1);
    polz_contriubutions=nan(numel(optical_freq),iimax);
    for ii=1:iimax
        qnum_excited_j_this=qnum_excited_j_vec(ii);
        freq_trans_this= freq_trans_vec(ii);
        first_part=( (-1)^qnum_excited_j_this  ) ...
            * Wigner6j(1,K,1,qnum_ground_j,qnum_excited_j_this,qnum_ground_j) ;
          
        polz_contriubutions(:,ii)= first_part ...
            * red_mat_elm_sq(ii) ...
            * (1/const.hb)* ( ...
                (1./(2*pi*freq_trans_this - 2*pi*optical_freq)) ...
                +(( (-1)^K )./(2*pi*freq_trans_this + 2*pi*optical_freq)) );
    end
    polz_k_array(:,K+1)= prefactor*sum(polz_contriubutions,2);
end

% equation 16

polz_scalar=polz_k_array(:,1) * 1/(sqrt(3*(2*qnum_ground_j+1)));
polz_vector=-polz_k_array(:,2) * sqrt(2*qnum_ground_j/(  (qnum_ground_j+1)*(2*qnum_ground_j+1)  ) );
polz_tensor=-polz_k_array(:,2) * sqrt(...
                                2*qnum_ground_j*(2*qnum_ground_j-1) ...
                                /( 3*(qnum_ground_j+1)*(2*qnum_ground_j+1)*(2*qnum_ground_j+3)  ) ...
                                );
if apply_polz_offset
    polz_scalar=polz_scalar+n_gt_3_polz*atom_u.polz;
end

higher_terms=[];
higher_terms.si_polz=[];
higher_terms.si_polz.scalar=polz_scalar;
higher_terms.si_polz.vector=polz_vector;
higher_terms.si_polz.tensor=polz_tensor;

end
