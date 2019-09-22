function str_out=num_with_sep(num_in,varargin)
% a decent implimentation if a little more complicated is also here
% https://au.mathworks.com/matlabcentral/fileexchange/52832-num2sepstr
% this implementation provides ability to handle the seperator after the decimal place


% examples
% num_with_sep(123456.012345)
% num_with_sep(123456.012345,'format','%.2f')
% num_with_sep(123456.012345,'separator',' ')
% num_with_sep(123456.012345,'add_sep_after_decimal',0)


%% parse the optional inputs
% format,add_sep_after_decimal,remove_trailing_zero,separator

is_c_logical=@(in) isequal(in,true) || isequal(in,false); %can x be cast as a logical
p = inputParser;
addParameter(p,'format','',@(x) ischar(x) || isstring(x) );
addParameter(p,'add_sep_after_decimal',true,is_c_logical);
addParameter(p,'remove_trailing_zero',nan,is_c_logical);
addParameter(p,'separator',',',@(x) (ischar(x) || isstring(x)) && numel(x)==1 );
parse(p,varargin{:});

format=p.Results.format;
add_sep_after_decimal=p.Results.add_sep_after_decimal;
remove_trailing_zero=p.Results.remove_trailing_zero;
num_sep_char=char(p.Results.separator);

if isempty(format)
    % stolen from % https://au.mathworks.com/matlabcentral/fileexchange/52832-num2sepstr
     if isinteger(num_in) || mod(round(num_in, 4), 1) == 0
        format = '%.0f';
    else
        format = '%.4f'; % 4 digits is the num2str default
     end
     default_format=true;
else
    default_format=false;
end

if isnan(remove_trailing_zero)
    if default_format
        remove_trailing_zero=true;
    else
        remove_trailing_zero=false;
    end
end

% solution from jan's comment here https://au.mathworks.com/matlabcentral/fileexchange/41520-insert-commas-into-numbers
num_str = sprintf(format, num_in); 
% we need to take out the core part with numbers
% eg '+123.456e-09' should return '123.456'
numeric_chars=isstrprop(num_str,'digit');
idx_dot=strfind(num_str,'.');
idx_dot=idx_dot(numeric_chars(idx_dot-1));
if ~isempty(idx_dot)
    idx_dot=idx_dot(1);
    numeric_chars(idx_dot)=true;
end
core_idx_start=find(numeric_chars,1, 'first');
core_idx_end=find(~numeric_chars(core_idx_start:end),1, 'first');
if isempty(core_idx_end)
    core_idx_end=numel(num_str);
else
    %subtract 2
    % one to account for core_idx_start starting at 1
    % one to account for the find being the first false value
    core_idx_end=core_idx_end+core_idx_start-2;
end


num_str_core= num_str(core_idx_start:core_idx_end);
num_str_xtras={num_str(1:core_idx_start-1),num_str(core_idx_end+1:end)};

decimal_place=strfind(num_str_core, '.');
% do not try to strip zeros from right if there is no decimal place 
if remove_trailing_zero && ~isempty(decimal_place) 
    num_str_core=strip(num_str_core,'right','0');
    decimal_place=strfind(num_str_core, '.');
end

if isempty(decimal_place)
    decimal_place=numel(num_str_core)+1;
end

num_str_core(2,(decimal_place-4):-3:1) = num_sep_char; 
if add_sep_after_decimal
    num_str_core(2,decimal_place+3:3:size(num_str_core,2)-1) = num_sep_char; 
end
num_str_core = transpose(num_str_core(num_str_core ~= char(0)));

str_out=[num_str_xtras{1},num_str_core,num_str_xtras{2}];

end