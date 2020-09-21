% what if only the large-R particles get pairs?

num_shots = 300;

mean_counts = 500;
count_var = 10;
num_counts = max(0,mean_counts + round(count_var*randn(num_shots,1)));

mean_wing_num = 300;
wing_num_var = 10;
num_wing_cts = max(mean_wing_num + round(wing_num_var*randn(num_shots,1)),0);
wing_size = 2;

mode = 'bimodal';
pair_frac = .01; %drop to .08 for realistic sim?
width = 0.1;
r_min = 0;


cli_header(1,'Synthesizing data...');
fake_counts = cell(1,num_shots);
for shot_idx = 1:num_shots
    
    % Generate a 3D normally distributed cloud 
    this_shot_counts = num_counts(shot_idx);
    init_counts = randn(this_shot_counts,3);
    
    if strcmp(mode,'unimodal')
        % Pick out a random subset to create pairs
        rand_idxs = randperm(this_shot_counts);
        num_pairs = round(pair_frac*this_shot_counts);
        % Make the back-to-back pairs with some noise
        buddy_particles = init_counts(rand_idxs(1:num_pairs),:);
        correlated_counts = - buddy_particles + width*randn(num_pairs,3);
        all_counts = [init_counts;correlated_counts];
        % set some unused values
    elseif strcmp(mode,'bimodal')        
        % Generate some extra counts from another distribution
        this_wing_num= num_wing_cts(shot_idx);
        wing_counts = wing_size*randn(this_wing_num,3);
        % select a random subset
        num_pairs = round(pair_frac*this_wing_num);
        rand_idxs = randperm(this_wing_num);
        % Make back-to-back pairs
        buddy_particles = wing_counts(rand_idxs(1:num_pairs),:);
        correlated_counts = - buddy_particles + width*randn(num_pairs,3);
        all_counts = [init_counts;wing_counts;correlated_counts];
    end

    fake_counts{shot_idx} = all_counts;
end
cli_header(2,'Done.');
% simulate the box
det_lims = [-1,1;-1,1;-5,5];
for ax = 1:3
   mask = all_counts(:,ax) > det_lims(ax,1) &  all_counts(:,ax) < det_lims(ax,2);
   all_counts = all_counts(mask);
end

stfig('Testing!');
clf;
subplot(2,2,1)
scatter3(all_counts(:,1),all_counts(:,2),all_counts(:,3));
xlim([-3,3])
ylim([-3,3])
zlim([-3,3])
for i=1:3
    subplot(2,2,i+1)
    histogram(all_counts(:,i),30);
end

%%
corr_opts.type='radial_bb';
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*5e-2; %only used for prewindow
rmax=5;
rmin=0;
corr_opts.redges=sqrt(linspace(rmin^2,rmax^2,600));
%corr_opts.redges=linspace(rmin,rmax,500);
corr_opts.rad_smoothing=0;
corr_opts.print_update = true;
corr_opts.plots = true;
corr_opts.do_pre_mask=false;
corr_opts.attenuate_counts=1;
corr_opts.one_d_smoothing=nan;
corr_opts.do_pre_mask=true;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=true;
corr_opts.low_mem=true;

% old code 27.137 s
tic
out=calc_any_g2_type(corr_opts,fake_counts);
toc