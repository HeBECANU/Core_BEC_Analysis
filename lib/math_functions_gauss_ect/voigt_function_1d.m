function v=voigt_function_1d(x,varargin)

% this function wraps up somf fuctions that calculate the voigt function
% into a slightly easier to use version where you input sigma and
% gamma of the gaussian and lorentzian respectively.
% voight calculators 
%       voigtf
%           https://au.mathworks.com/matlabcentral/fileexchange/50355-a-rational-approximation-of-the-voigt-function
%           in my basic testing this is able to produce a fractional error of about 2e-11
%       voigt_faddeeva
%           based on function in https://au.mathworks.com/matlabcentral/fileexchange/55334-voigt-model-fit
%           modified to use fadeva calculator https://au.mathworks.com/matlabcentral/fileexchange/47801-the-voigt-complex-error-function-second-version
%s
%
% mandatory inputs
%       x
% optional inputs 
%       sigma,gamma, mu,offset
% string pair args 
%       'norm' - can be 
%               'amp' - unit amplitude
%               'int' - unit integeral
%               'sum' - unit sum of elements
% 
% Bryce Henson 2020-04-22

x=col_vec(x);

p=inputParser;
valid_norm_type=@(x) any(strcmp(x,{'amp','int','sum'}));
valid_method_type=@(x) any(strcmp(x,{'approx','fadd','num'}));
addOptional(p,'sigma',1);
addOptional(p,'gamma',1);
addOptional(p,'mu',0.5);
addOptional(p,'amp',1);
addOptional(p,'offset',0);
addParameter(p,'norm','int',valid_norm_type);
addParameter(p,'method','approx',valid_method_type);
parse(p,varargin{:});

sigma = p.Results.sigma; 
gamma = p.Results.gamma; 
mu = p.Results.mu; 
amp = p.Results.amp; 
offset = p.Results.offset; 
norm_type=p.Results.norm;
method_type=p.Results.method;

% check that the sizes are either 
check_input_size(x,sigma)
check_input_size(x,gamma)
check_input_size(x,mu)
check_input_size(x,amp)
check_input_size(x,offset)


    omega_g=sigma*sqrt(2);

switch method_type
    case 'approx'
        % this aprox from for the voigt function
        % http://www.ccsenet.org/journal/index.php/jmr/article/view/62103
        % does not specity sigma and gamma seperately
        % so we must transform to the dimensionless form
        % see equation 6 of https://arxiv.org/pdf/0805.2274.pdf

        % which uses a different form of a gaussian so we define the gaussian omega from the gaussian sigma
    
        x_prime=(x-mu)./omega_g;
        y_prime=gamma/omega_g;
        v=(1/(omega_g*sqrt(pi)))*voigtf(x_prime,y_prime,2);
    case 'fadd'
        % basic idea from https://au.mathworks.com/matlabcentral/fileexchange/55334-voigt-model-fit
        % modified to use this fadeva calculator https://au.mathworks.com/matlabcentral/fileexchange/47801-the-voigt-complex-error-function-second-version
        % to see how this works check out http://www.ccsenet.org/journal/index.php/jmr/article/view/62103
        
        z = ((x-mu)+1i*gamma)/(omega_g);
        v = (1/(omega_g*sqrt(pi))) * real(fadf(z)); % Get Voigt from Faddeeva fn.
    case 'num'
        error('not yet implemented')
end

v=col_vec(v);

switch norm_type
    case 'amp'
        % this is rather crude but should work just fine
        % would nice to have some expression for the peak amplitude
        v=v/max(v(:));
    case 'sum'
        v=v/sum(v(:));
    case 'int'
         v=v;
end

c=v*amp +offset;



end


function check_input_size(x,var_to_check)
    if numel(var_to_check)~=1 && size(var_to_check)~=size(x)
        error('wrong size input must be the same size as x or scalar')
    end
end