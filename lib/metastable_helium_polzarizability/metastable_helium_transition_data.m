function transitions=metastable_helium_transition_data

global const

atom_u=[];
% atom_u.energy=(const.hb)^2/(const.me*(const.a0^2)); 
% atom_u.time=const.hb/atom_u.energy;
atom_u.energy=4.3597447222071e-18 ;%J

%% make a cell array of stucts with name, freq and osc strength
transitions=cell(0,0);



trans_23p0.name='2^{3}s_{1}-2^{3}p_{0}';
trans_23p0.type='E1'; % i think
trans_23p0.frequency_hz.val=276764094746.9e3;
trans_23p0.frequency_hz.unc=1.3e3;
trans_23p0.frequency_hz.ref.link='https://doi.org/10.1103/PhysRevLett.92.023001';
trans_23p0.frequency_hz.ref.note='table 1';
% trans_23p0.osc_strength.val=6.4175;
% trans_23p0.osc_strength.unc=nan;
% trans_23p0.osc_strength.ref.link='';
trans_23p0.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_23p0.osc_strength.ref.note='found from reduced matrix elements';
trans_23p0.osc_strength.units='in atomic energy units (hartree)';
reduced_matrix_elm=2.531782;
osc_strength=(trans_23p0.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_23p0.osc_strength.val=osc_strength;
trans_23p0.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_23p0);

trans_23p1.name='2^{3}s_{1}-2^{3}p_{1}';
trans_23p1.type='E1'; % i think
trans_23p1.frequency_hz.val=276764477805.0e3;
trans_23p1.frequency_hz.unc=0.9e3;
trans_23p1.frequency_hz.ref.link='https://doi.org/10.1103/PhysRevLett.92.023001';
trans_23p1.frequency_hz.ref.note='table 1';
% trans_23p1.osc_strength.val=19.253;
% trans_23p1.osc_strength.unc=nan;
% trans_23p1.osc_strength.ref.link='';
trans_23p1.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_23p1.osc_strength.ref.note='found from reduced matrix elements';
trans_23p1.osc_strength.units='in atomic energy units (hartree)';
reduced_matrix_elm=4.38477;
osc_strength=(trans_23p1.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_23p1.osc_strength.val=osc_strength;
trans_23p1.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_23p1);

trans_23p2.name='2^{3}s_{1}-2^{3}p_{2}';
trans_23p2.type='E1'; % i think
trans_23p2.frequency_hz.val=276732186818.4e3;
trans_23p2.frequency_hz.unc=1.5e3;
trans_23p2.frequency_hz.ref.link='https://doi.org/10.1103/PhysRevLett.92.023001';
trans_23p2.frequency_hz.ref.note='table 1';
% trans_23p2.osc_strength.val=32.0877;
% trans_23p2.osc_strength.unc=nan;
% trans_23p2.osc_strength.ref.link='http://dx.doi.org/10.1103/PhysRevA.88.052515';
% trans_23p2.osc_strength.ref.note='Table 1';
trans_23p2.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_23p2.osc_strength.ref.note='found from reduced matrix elements';
trans_23p2.osc_strength.units='atomic energy units (hartree)';
reduced_matrix_elm=5.66059;
osc_strength=(trans_23p2.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_23p2.osc_strength.val=osc_strength;
trans_23p2.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_23p2);


trans_33p0=[];
trans_33p0.name='2^{3}s_{1}-3^{3}p_{0}';
trans_33p0.type='E1'; % i think
trans_33p0.frequency_hz.val=770732839058e3;
trans_33p0.frequency_hz.unc=190e3;
trans_33p0.frequency_hz.ref.link='https://doi.org/10.1103/PhysRevLett.73.42';
% trans_33p0.osc_strength.val=0.2729;
% trans_33p0.osc_strength.unc=nan;
% trans_33p0.osc_strength.ref.link='http://dx.doi.org/10.1103/PhysRevA.88.052515';
% trans_33p0.osc_strength.ref.note='Table 1';
trans_33p0.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_33p0.osc_strength.ref.note='found from reduced matrix elements';
trans_33p0.osc_strength.units='in atomic energy units (hartree)';
reduced_matrix_elm=0.52478;
osc_strength=(trans_33p0.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_33p0.osc_strength.val=osc_strength;
trans_33p0.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_33p0);

trans_33p1=[];
trans_33p1.name='2^{3}s_{1}-3^{3}p_{1}';
trans_33p1.type='E1'; % i think
trans_33p1.frequency_hz.val=trans_33p0.frequency_hz.val - 8113.714e6;
trans_33p1.frequency_hz.unc=rssq([trans_33p0.frequency_hz.unc,0.028e6]);
trans_33p1.frequency_hz.ref.link='https://doi.org/10.1103/PhysRevLett.94.133001';
% trans_33p1.osc_strength.val=0.8188;
% trans_33p1.osc_strength.unc=nan;
% trans_33p1.osc_strength.ref.link='http://dx.doi.org/10.1103/PhysRevA.88.052515';
% trans_33p1.osc_strength.ref.note='Table 1';
trans_33p1.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_33p1.osc_strength.ref.note='found from reduced matrix elements';
trans_33p1.osc_strength.units='in atomic energy units (hartree)';
reduced_matrix_elm=0.90872;
osc_strength=(trans_33p1.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_33p1.osc_strength.val=osc_strength;
trans_33p1.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_33p1);


trans_33p2=[];
trans_33p2.name='2^{3}s_{1}-3^{3}p_{2}';
trans_33p2.type='E1'; % i think
trans_33p2.frequency_hz.val=trans_33p0.frequency_hz.val - 658.810e6;
trans_33p2.frequency_hz.unc=rssq([trans_33p1.frequency_hz.unc,0.018e6]);
trans_33p2.frequency_hz.ref.link='https://doi.org/10.1103/PhysRevLett.94.133001';
% trans_33p2.osc_strength.val=1.3646;
% trans_33p2.osc_strength.unc=nan;
% trans_33p2.osc_strength.ref.link='http://dx.doi.org/10.1103/PhysRevA.88.052515';
% trans_33p2.osc_strength.ref.note='Table 1';
trans_33p2.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_33p2.osc_strength.ref.note='found from reduced matrix elements';
trans_33p2.osc_strength.units='in atomic energy units (hartree)';
reduced_matrix_elm=1.17301;
osc_strength=(trans_33p2.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_33p2.osc_strength.val=osc_strength;
trans_33p2.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_33p2);


trans_43p0=[];
trans_43p0.name='2^{3}s_{1}-4^{3}p_{0}';
trans_43p0.type='E1'; % i think
trans_43p0.frequency_hz.val= const.c/318.86549229e-9;
trans_43p0.frequency_hz.unc=nan;
trans_43p0.frequency_hz.ref.link='https://www.researchgate.net/publication/246371868_High_Precision_Calculations_for_Helium';
trans_43p0.frequency_hz.ref.title='High precision calculations for helium, Drake';
trans_43p1.frequency_hz.ref.note='used values from nist spectrum database';
trans_43p0.frequency_hz.ref.type='theory';
% trans_43p0.osc_strength.val=1.3646;
% trans_43p0.osc_strength.unc=nan;
% trans_43p0.osc_strength.ref.link='http://dx.doi.org/10.1103/PhysRevA.88.052515';
% trans_43p0.osc_strength.ref.note='Table 1';
trans_43p0.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_43p0.osc_strength.ref.note='found from reduced matrix elements';
trans_43p0.osc_strength.units='in atomic energy units (hartree)';
reduced_matrix_elm=0.30039;
osc_strength=(trans_43p0.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_43p0.osc_strength.val=osc_strength;
trans_43p0.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_43p0);

trans_43p1=[];
trans_43p1.name='2^{3}s_{1}-4^{3}p_{1}';
trans_43p1.type='E1'; % i think
trans_43p1.frequency_hz.val= const.c/318.86661405e-9;
trans_43p1.frequency_hz.unc=nan;
trans_43p1.frequency_hz.ref.link='https://www.researchgate.net/publication/246371868_High_Precision_Calculations_for_Helium';
trans_43p1.frequency_hz.ref.title='High precision calculations for helium, Drake';
trans_43p1.frequency_hz.ref.note='used values from nist spectrum database';
trans_43p1.frequency_hz.ref.type='theory';
% trans_43p1.osc_strength.val=1.3646;
% trans_43p1.osc_strength.unc=nan;
% trans_43p1.osc_strength.ref.link='http://dx.doi.org/10.1103/PhysRevA.88.052515';
% trans_43p1.osc_strength.ref.note='Table 1';
trans_43p1.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_43p1.osc_strength.ref.note='found from reduced matrix elements';
trans_43p1.osc_strength.units='in atomic energy units (hartree)';
reduced_matrix_elm=0.52018;
osc_strength=(trans_43p1.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_43p1.osc_strength.val=osc_strength;
trans_43p1.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_43p1);

trans_43p2=[];
trans_43p2.name='2^{3}s_{1}-4^{3}p_{2}';
trans_43p2.type='E1'; % i think
trans_43p2.frequency_hz.val= const.c/318.86670551e-9;
trans_43p2.frequency_hz.unc=nan;
trans_43p2.frequency_hz.ref.link='https://www.researchgate.net/publication/246371868_High_Precision_Calculations_for_Helium';
trans_43p2.frequency_hz.ref.title='High precision calculations for helium, Drake';
trans_43p1.frequency_hz.ref.note='used values from nist spectrum database';
trans_43p2.frequency_hz.ref.type='theory';
% trans_43p2.osc_strength.val=1.3646;
% trans_43p2.osc_strength.unc=nan;
% trans_43p2.osc_strength.ref.link='http://dx.doi.org/10.1103/PhysRevA.88.052515';
% trans_43p2.osc_strength.ref.note='Table 1';
trans_43p2.osc_strength.ref.link='https://doi.org/10.1103/PhysRevA.93.052516';
trans_43p2.osc_strength.ref.note='found from reduced matrix elements';
trans_43p2.osc_strength.units='in atomic energy units (hartree)';
reduced_matrix_elm=0.67149;
osc_strength=(trans_43p2.frequency_hz.val*const.h/atom_u.energy) * (reduced_matrix_elm^2) * (2/(3*(2*1+1))) ;
trans_43p2.osc_strength.val=osc_strength;
trans_43p2.osc_strength.unc=nan;
transitions=cat(1,transitions,trans_43p2);


% Forbidden transitions


%might be better to use a combination of https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.111.013002 and the 1557 nm
%measurment instead of this direct measurment optical measurment for the frequency
% not sure where the reference for the oscillator strength has gone
% http://ffk2017.fuw.edu.pl/assets/1/3_Vassen.pptx but dead link
% will have to hunt through some theory papers for it
trans_21p1=[];
trans_21p1.name='2^{3}s_{1}-2^{1}p_{1}';
trans_21p1.type='E1+M2';
trans_21p1.einstein_a_coeff.val=1.548945;
trans_21p1.einstein_a_coeff.unc=nan;
trans_21p1.einstein_a_coeff.ref.link='https://arxiv.org/pdf/physics/0105110.pdf';
trans_21p1.frequency_hz.val=338133594.4e6;
trans_21p1.frequency_hz.unc=0.5e6;
trans_21p1.frequency_hz.ref.link='https://doi.org/10.1103/PhysRevLett.112.253002';
% now we need to turn this into an osc strength
% put all the ocnstants into a variable
rate_to_osc_strength_coef=(2*pi*(const.electron^2))/(const.epsilon0*const.me*(const.c^3));
multiplicity_ratio=1; %not sure if this is true
osc_strength=trans_21p1.einstein_a_coeff.val/(rate_to_osc_strength_coef*(trans_21p1.frequency_hz.val^2)*multiplicity_ratio);
trans_21p1.osc_strength.val=osc_strength;
trans_21p1.osc_strength.unc=nan;
trans_21p1.osc_strength.ref.link='derived from einstein_a_coeff';
transitions=cat(1,transitions,trans_21p1);

trans_21s0.name='2^{3}s_{1}-2^{1}S_{0}';
trans_21s0.type='M1'; %https://doi.org/10.1209%2Fepl%2Fi2006-10284-4
trans_21s0.einstein_a_coeff.val=9.1e-8;
trans_21s0.einstein_a_coeff.unc=nan;
trans_21s0.einstein_a_coeff.ref.link='https://doi.org/10.1209%2Fepl%2Fi2006-10284-4';
trans_21s0.einstein_a_coeff.ref.note='the paper "Helium 2 3S-2 1S metrology at 1.557 μm" quotes private coms with Pachucki ';
trans_21s0.frequency_hz.val=192510702145.6e3;
trans_21s0.frequency_hz.unc=1.8e3;
trans_21s0.frequency_hz.ref.link='https://doi.org/10.1126/science.1205163';
% now we need to turn this into an osc strength
% put all the ocnstants into a variable
rate_to_osc_strength_coef=(2*pi*(const.electron^2))/(const.epsilon0*const.me*(const.c^3));
multiplicity_ratio=1; %not sure if this is true
osc_strength=trans_21s0.einstein_a_coeff.val/(rate_to_osc_strength_coef*(trans_21s0.frequency_hz.val^2)*multiplicity_ratio);
trans_21s0.osc_strength.val=osc_strength;
trans_21s0.osc_strength.unc=nan;
trans_21s0.osc_strength.ref.link='derived from einstein_a_coeff';
transitions=cat(1,transitions,trans_21s0);


trans_33s1.name='2^{3}s_{1}-3^{3}S_{1}';
trans_33s1.type='M1'; %https://doi.org/10.1209%2Fepl%2Fi2006-10284-4
trans_33s1.einstein_a_coeff.val=6.484690e-9;
trans_33s1.einstein_a_coeff.unc=nan;
trans_33s1.einstein_a_coeff.ref.link='https://arxiv.org/abs/physics/0105110v1';
trans_33s1.einstein_a_coeff.ref.note='Table 1';
trans_33s1.frequency_hz.val=700939271e6;
trans_33s1.frequency_hz.unc=5e6;
trans_33s1.frequency_hz.ref.link='https://doi.org/10.1103/PhysRevLett.125.013002';
% now we need to turn this into an osc strength
% put all the ocnstants into a variable
rate_to_osc_strength_coef=(2*pi*(const.electron^2))/(const.epsilon0*const.me*(const.c^3));
multiplicity_ratio=1; %not sure if this is true
osc_strength=trans_33s1.einstein_a_coeff.val/(rate_to_osc_strength_coef*(trans_33s1.frequency_hz.val^2)*multiplicity_ratio);
trans_33s1.osc_strength.val=osc_strength;
trans_33s1.osc_strength.unc=nan;
trans_33s1.osc_strength.ref.link='derived from einstein_a_coeff';
transitions=cat(1,transitions,trans_33s1);



%%

% SetPrecision[{{, , 
% 
% 
% 
%     {f427, \
% ((0.1065403108*4.35974417*10^-18)/(6.62607004*10^-34)), "23s1-33s1"}},
%    50];
% 
% 
% 
% SetPrecision[{{0.2729, 770732839058*10^3, 
%     "23s1-33p0"}, {0.8188, (770732839058 - 811371)*10^3, 
%     "23s1-33p1"}, {1.3646, (770732839058 - 877252)*10^3, 
%     "23s1-33p2"}, {1.4878 10^-6, 338133594.4 10^6, 
%     "23s1-21p1"}, {6.4175, 276764094746.9*10^3, "23s1-23p0"}, {19.253,
%      276764477805.0*10^3, "23s1-23p1"}, {32.0877, 276732186818.4*10^3,
%      "23s1-23p2"}, {f1557, 192510702145.6*10^3, 
%     "23s1-21s0"}, {f427, \
% ((0.1065403108*4.35974417*10^-18)/(6.62607004*10^-34)), "23s1-33s1"}},
%    50];



% trans_template.name='';
% trans_template.osc_strength.val=nan;
% trans_template.osc_strength.unc=nan;
% trans_template.osc_strength.ref.link='';
% trans_template.osc_strength.ref.note='';
% trans_template.frequency_hz.val=nan;
% trans_template.frequency_hz.unc=nan;
% trans_template.frequency_hz.ref.link='';
% trans_template.osc_strength.ref.note='';
% transitions=cat(1,transitions,trans_template);


end