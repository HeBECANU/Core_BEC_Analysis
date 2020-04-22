function y=lorentzian_function_1d(x,varargin)
% mandatory inputs
%       x
% optional inputs 
%       gamma (half width at half maximum), mu (x center),offset
% string pair args 
%       'norm' - can be 
%               'amp' - unit amplitude
%               'int' - unit integeral
%               'sum' - unit sum of elements
% 
% Bryce Henson 2020-04-22

p=inputParser;
valid_norm_type=@(x) any(strcmp(x,{'amp','int','nint'}));
addOptional(p,'gamma',1);
addOptional(p,'mu',0.5);
addOptional(p,'amp',1);
addOptional(p,'offset',0);
addParameter(p,'norm','amp',valid_norm_type);
parse(p,varargin{:});
   
gamma = p.Results.gamma; 
mu = p.Results.mu; 
amp = p.Results.amp; 
offset = p.Results.offset; 
norm_type=p.Results.norm;

% check that the sizes are either 
check_input_size(x,gamma)
check_input_size(x,mu)
check_input_size(x,amp)
check_input_size(x,offset)


y = (gamma.^2)./((x-mu).^2 + gamma.^2);
    

switch norm_type
    case 'amp'
        y=y;
    case 'nint'
        y=y/trapz(x,y);
    case 'int'
         y=y/(gamma*pi);
end

y=y+offset;

end

function check_input_size(x,var_to_check)
    if numel(var_to_check)~=1 && size(var_to_check)~=size(x)
        error('wrong size input must be the same size as x or scalar')
    end
end