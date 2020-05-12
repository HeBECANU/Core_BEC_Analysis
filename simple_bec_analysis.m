%A tutorial on how to get up and running doing analysis on our BEC data
% lets begin
% it is common practice for the highest level file to have some user options in it
% for our uses this is a very reasonable practice, for bigger projects it could get a tad out of hand
%%%%------------------------start user var-------------------------------%
%first we need to point at the data directory
import_opts.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
%then what the MCP-TDC data files are called 
% the practice after 2016 is for it to just be 'd' eg d123.txt
import_opts.file_name='d';
%all the other feilds are optional, and have some sensible defaults
%in order of commonly used
%   the spatio-temporal limits on the data import
%import_opts.txylim=[[0,2];[-30e-3, 30e-3];[-30e-3, 30e-3]];     %tight XY lims to eliminate hot spot from destroying pulse widths
%   the rotation angle in radians
%import_opts.dld_xy_rot=0.61;
%import_opts.force_reimport=false;
%import_opts.force_forc=false;
%%%%------------------------done user var-------------------------------%

%% Setting Up The Enviorment

addpath('./lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path


% % find this .m file's path, this must be in the project root dir
% this_folder = fileparts(which(mfilename));
% % Add that folder plus all subfolders to the path.
% addpath(genpath(this_folder));%add all subfolders to the path to find genpath_exclude which should be in your project somewhere
% path_to_genpath=fileparts(which('genpath_exclude')); %save the dir that genpath_exclude is in
% path(pathdef) %clean up the path back to the default state to remove all the .git that were added
% addpath(this_folder)
% addpath(path_to_genpath)
% addpath(genpath_exclude(fullfile(this_folder,'lib'),'\.')) %add only a few directories,dont add hidden folders
% addpath(genpath_exclude(fullfile(this_folder,'dev'),'\.'))
% addpath(genpath_exclude(fullfile(this_folder,'bin'),'\.'))


hebec_constants %call the constants function that makes some globals


%% Importing Data
%now we can find all the data files in that directory
import_opts.shot_num=find_data_files(import_opts);
%and import them all into a structure, this function will also be cached in this directory so if you call the
%same import_opts it will just load the cache
data=import_mcp_tdc_data(import_opts);

%% Explore The Data Structure
% we have the fields
%            shot_num: [1×4140 double]
%            num_counts: [1×4140 double]
%            counts_txy: {1×4140 cell}
%            time_create_write: [4140×2 double]
% the one we care the most about is counts_txy
% we can look at the first file's distribution of detections in time

stfig(1); % quiet figure that wont steal focus, 
          % also has the option to show what function called it which is useful when you have lots of functions making figs
clf
% lets make a basic histogram of the time arrivals of atoms
histogram(data.counts_txy{1}(:,1),linspace(0,1,1e3))
xlabel('time(s)')
ylabel('counts')

%% Combining multiple shots
% now lets plot the same thing but for the first 10 shots combined
txy_multiple_shots=vertcat(data.counts_txy{1:3});
% lets do a better histogram, if you take a histogram with bins that are smaller than normal and then smooth(convolve
% with a guassian) you get something called a kernel density estimation which reduced the problems histograms have with
% hiding features (eg a big dip in counts that is smaller than the bin width)
% i have made a function that does this
% set up the options for the smooth histograming function
out_struct=smooth_hist(txy_multiple_shots(:,1),'lims',[0,2],'sigma',1e-4,'bin_factor',10);
stfig('count rate TOF & spectrum','add_stack',1); %this time we will give the figure a name and prepend the function that called it
clf
subplot(2,1,1)
plot(out_struct.bin.centers,out_struct.count_rate.smooth*1e-3)
ylabel('count rate (khz)')
xlabel('time (s)')
xlim([0.4,0.7])
subplot(2,1,2)

%% Investigating the spectrum
% use a fft to find what frequency components are
% this fft_tx function makes taking a fft much easier, it handles the padding and windowing which can be annoying to get
% right. It spits out the bin frequencies and amplitudes.
fft_out=fft_tx(out_struct.bin.centers,out_struct.count_rate.smooth,'window','chebyshev','win_param',{300},'padding',10);

subplot(2,1,2)
cla
semilogy(fft_out(1,:),abs(fft_out(2,:)))
xlim([0,1000])
%ylim([1e2,1e6])

%% Finding the main components in the spectrum
% this function will find the dominant freqency componets in a timeseries
dom_opt=[];
dom_opt.num_components=10;
[components,details]=dominant_freq_components(out_struct.bin.centers,out_struct.count_rate.smooth,dom_opt);

%semilogy(details.fft_dat(1,:),abs(details.fft_dat(2,:)))

% That should give you a decent start
% for more advanced applications you will want to loop over each shot
% have a look through the lib folder for more things
