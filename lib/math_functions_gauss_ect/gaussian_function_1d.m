function y=gaussian_function_1d(x,varargin)
% mandatory inputs
%       x
% optional inputs 
%       sigma, mu,amp,offset
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
addParameter(p,'derivative',0,@(x) any(x==[0,1,2,3,4,5]))
parse(p,varargin{:});
   
sigma = p.Results.sigma; 
mu = p.Results.mu; 
amp = p.Results.amp; 
offset = p.Results.offset; 
norm_type=p.Results.norm;
derivative_order=p.Results.derivative;

% check that the sizes are either 
check_input_size(x,sigma)
check_input_size(x,mu)
check_input_size(x,amp)
check_input_size(x,offset)


    
switch derivative_order
    % for derivations see mathematica document
    case 0
        y = exp(-(1/2)*((x-mu)./sigma).^2);
    case 1
        y = -( (x-mu) ./ (sigma.^2) ).*exp(-(1/2)*((x-mu)./sigma).^2);
    case 2
    	xdiff=x-mu;  
        exp_part=exp(-(1/2)*((xdiff)./sigma).^2);
        % not sure if this is the best way to get numerical stability
        y= +exp_part .* ( (xdiff.^2)./(sigma.^4) )   ...
           -exp_part .* (1./( sigma.^2) ) ;
    case 3
        xdiff=x-mu;  
        exp_part=exp(-(1/2)*((xdiff)./sigma).^2);
        y=-exp_part .* ( (xdiff.^3)./(sigma.^6) ) ...
          +exp_part .* (3*(xdiff )./( sigma.^4 ) )  ;
    case 4
        xdiff=x-mu;  
        exp_part=exp(-(1/2)*((xdiff)./sigma).^2);
        y=+exp_part .* (( xdiff.^4 )./sigma.^8) ... 
          -exp_part .* ((6*xdiff.^2)./(sigma.^6) ) ...
          +exp_part .* (3./(sigma.^4) )  ;
    case 5
        xdiff=x-mu;  
        exp_part=exp(-(1/2)*((xdiff)./sigma).^2);
        y=-exp_part .* ( (xdiff.^5)./(sigma.^10) ) ...
          +exp_part .* 10.*( (xdiff.^3)./(sigma.^8) ) ...
          -exp_part.* 15.* (xdiff./(sigma.^6) );
end


switch norm_type
    case 'amp'
        y=y;
    case 'nint'
        if derivative_order~=0
            warning('numeric integeration normalization with the derviate option may not be what you want')
        end
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