function handle_out = stfig(fig_name_num_obj,varargin)
% stfig - silent text figure
% open and select figures based on text names withoug annoying focus-theft.
% built on sfig by Daniel Eaton, 2005
%
% Syntax:   fig_obj=stfig(in)
%           fig_obj=stfig('my new fig name')
%           fig_obj=stfig('my new fig name','add_stack',1) %add the function stack to the name
%           fig_obj=stfig('my new fig name','show_handle',1) %show the numeric fig handle eg 'Figure 5: my new fig name'.
%           fig_obj=stfig('my new fig name','return_cwin',1) %show the numeric fig handle eg 'Figure 5: my new fig name'.
%
%
%           stfig % create a new figure with the stack trace as the name (no numeric fig name)
%           stfig([]) %create a new numeric figure, with ('show_handle'=true) equivelent to figure
%           stfig([],'return_cwin',1) %create a new numeric figure with ('show_handle'=true) and return focus to the
%                                          command line
%           stfig('','add_stack',1) equivelent to stfig
%
% Inputs:
%   fig_name_num_obj - [fig text,fig number,fig object]
%                    -fig text, char or string input of the figure name that exists or you want to create eg, 'my new fig'
%                    -fig number, figure number that exists or you want to create
%                    -fig object, a valid fig object
%                    -empty or NaN, create a new numeric figure
%   Optional Name Value Pairs:
%       'add_stack' - logical[default=false], add the function stack to the figure name eg 
%                                      test_stfig_text: test_stfig_mock_call_fun: my fig name
%       'show_handle'       - logical[default=false], show the figure number eg 'Figure 5: my new fig name'.
%       'return_cwin'    - logical[default=false], move focus back to the comand window after figure is created
% Outputs:
%    fig_obj - an object of the class matlab.ui.Figure
%
% Example: 

%         
% Other m-files required: none
% Also See: test_stfig_text,sfigure(https://au.mathworks.com/matlabcentral/fileexchange/8919-smart-silent-figure)
% Subfunctions: make_fun_stack_str,clip_str
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%    - acessing figures with same label after function stack
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com %YOU MUST INCLUDE '[matlab][stfig]' in the subject line OR I WILL NOT REPLY
% Last revision:2019-05-16


%------------- BEGIN USER VAR ------------
max_text_width=70;
fun_stack_msb_left=1; %should the stack text have the most signifigant function to the left(true)
add_str_left=0; %should the input text be added to the left(true) or right (false) Sensible to be ~fun_stack_msb_order
clip_str_left=1; %clip the end(right)(true) or start(left)(false) of the title text. Sensible to be ~add_str_left
sep_str=': ';

%------------- END USER VAR --------------
%------------- BEGIN CODE ----------------

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%parse user inputs
p = inputParser;
is_c_logical=@(in) isequal(in,true) || isequal(in,false); %can x be cast as a logical
addParameter(p,'add_stack',false,is_c_logical);
addParameter(p,'show_handle',false,is_c_logical);
addParameter(p,'return_cwin',false,is_c_logical);
parse(p,varargin{:});

add_stack_to_text=p.Results.add_stack;
show_handle=p.Results.show_handle;
focus_back_to_cwin=p.Results.return_cwin;
% use logical shortcuting to prevent execution of isnan(fig_name_num_obj) if fig_name_num_obj is a char vec
% 
if nargin>=1 
    if ~isempty(fig_name_num_obj)
        if ischar(fig_name_num_obj)||isstring(fig_name_num_obj)
            if nargin<3 || isempty(show_handle)
                show_handle=false;
            end
            if nargin<2 || isempty(add_stack_to_text)
                add_stack_to_text=false;
            end
            if add_stack_to_text
                fun_stack_str=make_fun_stack_str(fun_stack_msb_left);
                if add_str_left
                    fig_str=[fig_name_num_obj,sep_str,fun_stack_str];
                else
                    fig_str=[fun_stack_str,sep_str,fig_name_num_obj];
                end
            else
                fig_str=fig_name_num_obj;
            end
            fig_str=clip_str(fig_str,max_text_width,clip_str_left);

            % TODO: save the name that the figure was called with to the user data so that it can be found across 
            % different function stacks, fun_stack_str . Output all the figures and extract out the UserData
            %find all the existing figures
            %exisiting_figs=findall(groot, 'Type', 'figure');  
            %strings=arrayfun(@(x) x.Name,exisiting_figs,'UniformOutput',false);
            %match_idx=find(strcmp(strings,fig_str));

            exisiting_figs=findall(groot, 'Type', 'figure','Name',fig_str);
            if numel(exisiting_figs)>1
                warning('multiple figs with that name exist, will use the first instance')
            end  
            if numel(exisiting_figs)==1 %does the figure exist?
                set(0, 'CurrentFigure', exisiting_figs(1).Number)
                %TODO: this will change if the default value for show_handle is used, should check if the optional input
                %was passed
                exisiting_figs(1).NumberTitle=show_handle; 
                handle_out=exisiting_figs(1);
            else
                if show_handle
                    handle_out=figure('Name',fig_str);
                else
                    handle_out=figure('Name',fig_str,'NumberTitle','off');
                end
                set(handle_out,'color','w')
%                 set(handle_out,'DefaultAxesXLabelFontSize','20')
                if focus_back_to_cwin
                    drawnow
                    commandwindow
                end
            end
        elseif isnumeric(fig_name_num_obj)
            if round(fig_name_num_obj)~=fig_name_num_obj
                error('input must be integer')
            end
            if ishandle(fig_name_num_obj)
                set(0, 'CurrentFigure', fig_name_num_obj);
                handle_out=gcf;
            else
                handle_out = figure(fig_name_num_obj);
                set(handle_out,'color','w');
                if focus_back_to_cwin
                    drawnow
                    commandwindow
                end
            end
        elseif strcmp(class(fig_name_num_obj),'matlab.ui.Figure')
            if isvalid(fig_name_num_obj)
                set(0, 'CurrentFigure', fig_name_num_obj.Number)
                handle_out=fig_name_num_obj;
            else
                warning('fig object invalid, fig may have been closed, calling making new figure')
                handle_out=stfig([], varargin{:});
            end
        else
            warning('invalid input as first argument must be string,int, or figure obj')
        end
    else %isempty(fig_name_num_obj) case
        handle_out = figure;
        set(handle_out,'color','w');
        if focus_back_to_cwin
            drawnow
            commandwindow
        end
    end
else %no input arguments case
	%if no arguments call the function but set the fig name to be the stack trace
    %need to do make_fun_stack_str here or if we did stfig('','add_stack',1) then the stack would be one deeper 
    fun_stack_str=make_fun_stack_str(fun_stack_msb_left);
    handle_out=stfig(fun_stack_str);
end
end


function str_out=make_fun_stack_str(fun_stack_msb_left)
% fun_stack_msb_left %should the stack text have the most signifigant function to the left(true) or right(false)

sep_str=': ';

if nargin==0
    fun_stack_msb_left=1; 
end

fun_stack_obj=dbstack;
if fun_stack_msb_left
    fun_stack_obj=flipud(fun_stack_obj);
    fun_stack_txt=sprintf(['%s',sep_str],fun_stack_obj(1:end-2).name);
else
    fun_stack_txt=sprintf(['%s',sep_str],fun_stack_obj(3:end).name);
end
str_out=fun_stack_txt(1:end-numel(sep_str));   
end


function str_out=clip_str(str_in,max_text_width,clip_left)
% clip_right clip the right side of the text(true) or left(false)

cont_str='...';
str_in_len=numel(str_in);
if str_in_len>max_text_width
    clip_str_len=max_text_width-numel(cont_str);
    if clip_left
    	str_out=str_in(end-clip_str_len:end);
    	str_out=['...',str_out];
    else
        str_out=str_in(1:clip_str_len);
        str_out=[str_out,'...'];
    end
else
    str_out=str_in;
end

end