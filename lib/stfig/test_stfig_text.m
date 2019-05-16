%test stfig
% walk through use cases
close('all')

% open a figure with a text name
stfig('a noisy sine wave');
plot(sin(linspace(0,1,1e3))+rand(1,1e3))

% open another figure with the function stack appended to the name
stfig('another noisy plot','add_stack',1);
plot(rand(1,1e3))

test_stfig_mock_call_fun('my fig from a function');

%update a previously opened figure using the name
pause(1)
stfig('a noisy sine wave');
plot(-sin(linspace(0,1,1e3))+rand(1,1e3));

% make the handle name visible
pause(1)
fig_obj=stfig('a noisy sine wave','show_handle',1);
plot(-sin(linspace(0,1,1e3))+rand(1,1e3))

%% uisng handles
%plot something
fig_obj=stfig('a noisy sine wave');
plot(-sin(linspace(0,1,1e3))+rand(1,1e3))

% do something else
stfig('another noisy plot');
plot(rand(1,1e3))
% acess the figure using the first figure object that was output
stfig(fig_obj);
pause(0.5)
plot(sin(linspace(0,1,1e3))+rand(1,1e3))

% gracefully handles being passed a closed figure
pause(0.5)
close(fig_obj)
stfig(fig_obj);


%% ================limitations====================
%% cant find figures across add_stack_to_text options
% a figure opened with the add_stack_to_text=true
stfig('another noisy plot',1)
plot(rand(1,1e3))
% has a different name than if you call it with add_stack_to_text=false
stfig('another noisy plot',0)
plot(rand(1,1e3))
% this makes acessing figures across functions near impossible so if this is your use case you should aviod this use.
% future additions could use the UserData property of the figure to store the plain name to allow it to be
% found

%% no direct mechanism to change figure name
% however matlab has a very acessible way to do this with the figure object that
% the stfig funtion returns
fig_obj=stfig('a noisy sine wave');
fig_obj.Name='changed name';

%% speed
% searching through all the figures takse matlab >15ms
tic
fig_obj=stfig('a noisy sine wave');
toc
% however once opened using the fig_obj to set the current figure is <1ms
tic
stfig(fig_obj);
toc
% if you dont care about dunamicaly naming things then you can use the figure number directly, this against the spirit
% of this function but can be handy in dev
tic
stfig(3);
toc

