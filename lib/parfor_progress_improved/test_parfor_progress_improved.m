% basic usage

% only update some of the time so it wont slow down some fast inner loop
tic
updates=100;
itt = 100000;
update_frac=itt/updates;
parfor_progress_imp(updates);
parfor i=1:itt
    pause(10*rand/itt); % Replace with real code
    if rand<(1/update_frac)
        parfor_progress_imp;
    end
end
parfor_progress_imp(0);
toc


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
toc

%%

tic 
itt = 100;
parfor_progress(itt);
parfor i=1:itt
    parfor_progress;
end
parfor_progress_imp(0);
toc


