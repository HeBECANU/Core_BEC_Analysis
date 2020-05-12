% basic usage

% only update some of the time so it wont slow down some fast inner loop
tic
updates=100;
itt = 100000;
update_frac=updates/itt;
parfor_progress_imp(updates);
num_updates=0;
pause(1)
for i=1:itt
    pause(0.1*rand/itt); % Replace with real code
    if rand<(update_frac*1.3)
        parfor_progress_imp;
        num_updates=num_updates+1;
        pause(0.0002)
        if num_updates>0.95*updates
            pause(0.1)
        end
    end
end
parfor_progress_imp(0);
toc
num_updates



%%
% can also be used in nested loops

tic
updates=100;
itt_outer = 100000;
itt_inner=10;
update_frac=(itt*itt_inner)/updates;
parfor_progress_imp(updates);
parfor i=1:itt_outer
    for x=1:itt_inner
        pause(10*rand/itt); % Replace with real code
        if rand<(1/update_frac)
            parfor_progress_imp;
        end
    end
end
parfor_progress_imp(0);
toc

%% and you can deterministicly update instead of randomly
tic
updates=100;
itt_outer = 100;
itt_inner=100;
update_frac=(itt*itt_inner)/updates;
parfor_progress_imp(itt_outer*itt_inner/30);
for i=1:itt_outer
    for x=1:itt_inner
        pause(0.001); % Replace with real code
        if mod(x,30)==0
            parfor_progress_imp;
        end
    end
end
parfor_progress_imp(0);
toc

%% Benchamrk against original
tic 
itt = 100;
parfor_progress_imp(itt);
parfor i=1:itt
    parfor_progress_imp;
end
parfor_progress_imp(0);
time_imp=toc;
time_per_update_imp=time_imp/itt;
fprintf('time per update %.2e\n',time_per_update_imp)
%%

tic 
itt = 100;
parfor_progress(itt);
parfor i=1:itt
    parfor_progress;
end
parfor_progress_imp(0);
time_old=toc;
time_per_update_old=time_old/itt;
fprintf('time per update %.2e\n',time_per_update_old)


