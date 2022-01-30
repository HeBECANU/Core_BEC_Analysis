example_fit_pal_atom_number

% as ripped from tune out code


%% FITTING THE ATOM NUMBER
%use the inital few atom laser pulses in order to determine the atom number
%not of that much benifit TBH
anal_opts.atom_num_fit=[];
anal_opts.atom_num_fit.pulses=[2,20]; %min,max index of pulses
anal_opts.atom_num_fit.plot.each_shot=false;
anal_opts.atom_num_fit.plot.history=false;
anal_opts.atom_num_fit.qe=anal_opts.global.qe;

data.num_fit=fit_atom_number(anal_opts.atom_num_fit,data);
