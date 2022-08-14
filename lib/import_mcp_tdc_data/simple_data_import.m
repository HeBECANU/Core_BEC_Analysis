%% simple data import
opts.import.dir = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20220628_thermal_atom_laser';   %the directory where you are keeping the data
opts.import.cache_save_dir = fullfile(opts.import.dir, 'cache', 'import\');
opts.import.force_reimport = true;  %do you want to reimport the data or use a cache
opts.import.force_cache_load = ~opts.import.force_reimport;
%% Import parameters
tmp_xlim=[-35e-3, 35e-3];   %XY lims of imported data
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,6]; %time limits of imported data
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

        
%% actually import the data        
[data, ~] = import_mcp_tdc_data(opts.import);

%% mask txy
%here we do some further masking to select only one section of the data
shot = 5;% shot to consider
txylim = [3.1,3.2;
    -35e-3,35e-3;
    -35e-3,35e-3];
masked_data = masktxy_square(data.counts_txy{shot},txylim);
%% plot some data
figure(11)
histogram(masked_data(:,1),100)
xlabel('arrival time')
ylabel('number of counts')