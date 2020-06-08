function y=lorentzian_function_1d(x,varargin)
% mandatory inputs
%       x
% optional inputs 
%       gamma (half width at half maximum), mu (x center),amp,offset
% string pair args 
%       'norm' - can be 
%               'amp' - unit amplitude
%               'int' - unit integeral
%               'sum' - unit sum of elements
% 
% lor_fun=@(b,x) b(3)*(b(1).^2)./((x-b(2)).^2 + b(1).^2)+b(4);
% cof_names={'gamma','mu','amp','offset'};
% Bryce Henson 2020-04-22

p=inputParser;
valid_norm_type=@(x) any(strcmp(x,{'amp','int','nint'}));
addOptional(p,'gamma',1);
addOptional(p,'mu',0.5);
addOptional(p,'amp',1);
addOptional(p,'offset',0);
addParameter(p,'norm','amp',valid_norm_type);
addParameter(p,'derivative',0,@(x) any(x==[0,1,2]))
parse(p,varargin{:});
   
gamma = p.Results.gamma; 
mu = p.Results.mu; 
amp = p.Results.amp; 
offset = p.Results.offset; 
norm_type=p.Results.norm;
derivative_order=p.Results.derivative;

% check that the sizes are either 
check_input_size(x,gamma)
check_input_size(x,mu)
check_input_size(x,amp)
check_input_size(x,offset)

switch derivative_order
    case 0
        y = (gamma.^2)./((x-mu).^2 + gamma.^2);
    case 1
        y = -(2*(x-mu).*gamma.^2)./(  ((x-mu).^2 + gamma.^2).^2  );
    case 2
      xdiff=x-mu;  
      y=gamma.^2 * (   8*(xdiff.^2)./( (xdiff.^2+gamma.^2).^3  )   - ...
                       2./( (xdiff.^2+gamma.^2).^2  )  );
end

switch norm_type
    case 'amp'
        y=y;
    case 'nint'
        y=y/trapz(x,y);
    case 'int'
         y=y/(gamma*pi);
end
y=y*amp+offset;

end

function check_input_size(x,var_to_check)
    if numel(var_to_check)~=1 && size(var_to_check)~=size(x)
        error('wrong size input must be the same size as x or scalar')
    end
end