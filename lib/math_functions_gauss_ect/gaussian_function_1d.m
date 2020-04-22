function y=gaussian_function_1d(x,varargin)
% mandatory inputs
%       x
% optional inputs 
%       sigma, mu,offset
% string pair args 
%       'norm' - can be 
%               'amp' - unit amplitude
%               'int' - unit integeral
%               'sum' - unit sum of elements
% 
% Bryce Henson 2020-04-22


p=inputParser;
valid_norm_type=@(x) any(strcmp(x,{'amp','int','nint'}));
addOptional(p,'sigma',1);
addOptional(p,'mu',0.5);
addOptional(p,'amp',1);
addOptional(p,'offset',0);
addParameter(p,'norm','amp',valid_norm_type);
parse(p,varargin{:});
   
sigma = p.Results.sigma; 
mu = p.Results.mu; 
amp = p.Results.amp; 
offset = p.Results.offset; 
norm_type=p.Results.norm;

% check that the sizes are either 
check_input_size(x,sigma)
check_input_size(x,mu)
check_input_size(x,amp)
check_input_size(x,offset)


    
y =exp(-(1/2)*((x-mu)./sigma).^2);

switch norm_type
    case 'amp'
        y=y;
    case 'nint'
        y=y/trapz(x,y);
    case 'int'
         y=y/(sigma*sqrt(2*pi));
end

y=y.*amp+offset;

end

function check_input_size(x,var_to_check)
    if numel(var_to_check)~=1 && size(var_to_check)~=size(x)
        error('wrong size input must be the same size as x or scalar')
    end
end