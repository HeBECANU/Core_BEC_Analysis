% test sctipt for min_diff
% TODO
% 
% Other m-files required: isalmost https://au.mathworks.com/matlabcentral/fileexchange/15816-isalmost


%% Test basic usage

fprintf('======testing basic usage \n')

dat_elements=5;
min_diff_test_dat=1e-6;
logic_str = {'FAIL', 'pass'};
data_in=linspace(0,1,dat_elements-1);
%test with column vector
data_in=[data_in,1+min_diff_test_dat]';
%randomly permute the order
rand_order=randperm(numel(data_in));
data_in=data_in(rand_order);
eval_timer=tic;
[min_diff_val,min_diff_idxs]=min_diff_vec(data_in');
eval_time=toc(eval_timer);
%find where the last and second last points were put in the ransomized data list
idx_ans=sort([find(rand_order==dat_elements-1),find(rand_order==dat_elements)]);
idx_sort=sort(min_diff_idxs); %sort the output of min_diff for compare purposes
fprintf('TEST: index correct: %s\n',logic_str{1+isequal(idx_ans,idx_sort)})
fprintf('TEST: diff correct : %s\n',logic_str{1+isalmost(min_diff_test_dat,min_diff_val,eps)})
fprintf('INFO: eval time: %.3f ms\n',eval_time*1e3)


%% test tolerance

fprintf('======testing diff tolerance function \n')

dat_elements=5;
min_diff_test_dat=1e-6;
diff_tol=1e-9;
data_in=linspace(0,1,dat_elements-1);
%test with column vector
data_in=[data_in,1+min_diff_test_dat]; %add in the difference that we want to find
%add in some elements with differences smaller than tolerance
data_in=[data_in,-0.5*diff_tol];
%randomly permute the order
rand_order=randperm(numel(data_in));
data_in_rand=data_in(rand_order);
eval_timer=tic;
[min_diff_val,min_diff_idxs]=min_diff_vec(data_in_rand',diff_tol);
eval_time=toc(eval_timer);
%find where the last and second last points were put in the ransomized data list
idx_ans=sort([find(rand_order==dat_elements-1),find(rand_order==dat_elements)]);
idx_ans=idx_ans(:); %deal with things being row or col vec 
idx_sort=sort(min_diff_idxs(:)); %sort the output of min_diff for compare purposes
fprintf('TEST: index correct: %s\n',logic_str{1+isequal(idx_ans,idx_sort)})
fprintf('TEST: diff correct : %s\n',logic_str{1+isalmost(min_diff_test_dat,min_diff_val,eps)})
fprintf('INFO: eval time    : %.3f ms\n',eval_time*1e3)




%% 
fprintf('======comparing speed with brute algo \n')
fprintf('======Brute algo\n')
dat_elements=1e4;
min_diff_test_dat=1e-6;
data_in=linspace(0,1,dat_elements-1);
%test with column vector
data_in=[data_in,1+min_diff_test_dat]';
%randomly permute the order
rand_order=randperm(numel(data_in));
data_in=data_in(rand_order);
eval_timer=tic;
[min_diff_val,min_diff_idxs]=min_diff_vec_brute(data_in');
eval_time_brute=toc(eval_timer);
%find where the last and second last points were put in the ransomized data list
idx_ans=sort([find(rand_order==dat_elements-1),find(rand_order==dat_elements)]);
idx_sort=sort(min_diff_idxs); %sort the output of min_diff for compare purposes
fprintf('TEST: index correct: %s\n',logic_str{1+isequal(idx_ans,idx_sort)})
fprintf('TEST: diff correct : %s\n',logic_str{1+isalmost(min_diff_test_dat,min_diff_val,eps)})
fprintf('INFO: eval time    : %.3f ms\n',eval_time_brute*1e3)
fprintf('======Search algo \n')
eval_timer=tic;
[min_diff_val,min_diff_idxs]=min_diff_vec(data_in');
eval_time_search=toc(eval_timer);
%find where the last and second last points were put in the ransomized data list
idx_ans=sort([find(rand_order==dat_elements-1),find(rand_order==dat_elements)]);
idx_sort=sort(min_diff_idxs); %sort the output of min_diff for compare purposes
fprintf('TEST: index correct: %s\n',logic_str{1+isequal(idx_ans,idx_sort)})
fprintf('TEST: diff correct : %s\n',logic_str{1+isalmost(min_diff_test_dat,min_diff_val,eps)})
fprintf('INFO: eval time    : %.3f ms\n',eval_time_search*1e3)
fprintf('INFO: rel eval time: %.1f\n',eval_time_brute/eval_time_search)

%% compare method scaling
%%
%Should find the time for counting and returning as seprate plots
fprintf('======comparing scaling with brute algo \n')


points_list=round(logspace(0.5,4,400)); %values of n to investigate
points_list=unique(points_list);
min_diff_test_dat=1e-6;
%set up the output arrays
iimax=numel(points_list);
time_sort_method=nan(iimax);
time_brute_method=nan(iimax);

last_update=posixtime(datetime('now')); %time for updating plots every few seconds

tic
sfigure(1);
clf
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
plot_colors=parula(5+1); %padd the colors to avoid yellow

fprintf('  \n%03u',0) 
for ii=1:iimax
fprintf('\b\b\b%03u',ii)  


dat_elements=points_list(ii);
data_in=linspace(0,1,dat_elements-1);
%test with column vector
data_in=[data_in,1+min_diff_test_dat]';
%randomly permute the order
rand_order=randperm(numel(data_in));
data_in=data_in(rand_order);

eval_timer=tic;
[min_diff_val,min_diff_idxs]=min_diff_vec_brute(data_in');
time_brute_method(ii)=toc(eval_timer);
%find where the last and second last points were put in the ransomized data list
idx_ans=sort([find(rand_order==dat_elements-1),find(rand_order==dat_elements)]);
if ~(isequal(idx_ans,min_diff_idxs)) || ~isalmost(min_diff_test_dat,min_diff_val,eps)
    error('not equal')
end


eval_timer=tic;
[min_diff_val,min_diff_idxs]=min_diff_vec(data_in');
time_sort_method(ii)=toc(eval_timer);
%find where the last and second last points were put in the ransomized data list
idx_ans=sort([find(rand_order==dat_elements-1),find(rand_order==dat_elements)]);
if ~(isequal(idx_ans,min_diff_idxs)) || ~isalmost(min_diff_test_dat,min_diff_val,eps)
    error('not equal')
end


ptime=posixtime(datetime('now'));
if ptime-last_update>2 || ii==iimax
    sfigure(1);
    loglog(points_list,...
        time_brute_method,'color',plot_colors(1,:))
    hold on
    loglog(points_list,...
        time_sort_method,'color',plot_colors(3,:))
    legend('brute','unsorted',...
        'Location','northwest')
    hold off
    xlabel('Vector Size(n)');
    ylabel('Execution Time');
    title('Return smallest diff(and idx) between elements')

    pause(1e-6)
    last_update=ptime;

end

end
figure(1)
fprintf('\n') 

saveas(gcf,'fig1.png')


