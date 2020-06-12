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
addParameter(p,'derivative',0,@(x) any(x==[0,1,2]));
parse(p,varargin{:});

sigma = p.Results.sigma; 
gamma = p.Results.gamma; 
mu = p.Results.mu; 
amp = p.Results.amp; 
offset = p.Results.offset; 
norm_type=p.Results.norm;
method_type=p.Results.method;
derivative_order=p.Results.derivative;

% check that the sizes are either 
check_input_size(x,sigma)
check_input_size(x,gamma)
check_input_size(x,mu)
check_input_size(x,amp)
check_input_size(x,offset)

sigma=abs(sigma);
gamma=abs(gamma);


% if abs(sigma/gamma)<1e-3 && abs(sigma/gamma)>1e-6
%     method_type='fadd';
%     warning('sigma is too small using fadd method')
% end
switch derivative_order
    case 0
        if abs(sigma/gamma)<1e-6
           warning('using lorentzian because sigma is so small') 
           v=lorentzian_function_1d(x,gamma,mu,amp,0,'norm','amp');
        elseif abs(gamma/sigma)<1e-6
            warning('using gaussian because gamma so small')
            v=gaussian_function_1d(x,sigma,mu,amp,0,'norm','amp');
        else
            switch method_type
                case 'approx'
                    % this aprox from for the voigt function
                    % http://www.ccsenet.org/journal/index.php/jmr/article/view/62103
                    % does not specity sigma and gamma seperately
                    % so we must transform to the dimensionless form
                    % see equation 6 of https://arxiv.org/pdf/0805.2274.pdf

                    % which uses a different form of a gaussian so we define the gaussian omega from the gaussian sigma
                    omega_g=sigma*sqrt(2);
                    x_prime=(x-mu)./omega_g;
                    y_prime=gamma/omega_g;
                    v=(1/(omega_g*sqrt(pi)))*voigtf(x_prime,y_prime,2);
                case 'fadd'
                    % basic idea from https://au.mathworks.com/matlabcentral/fileexchange/55334-voigt-model-fit
                    % modified to use this fadeva calculator https://au.mathworks.com/matlabcentral/fileexchange/47801-the-voigt-complex-error-function-second-version
                    % to see how this works check out http://www.ccsenet.org/journal/index.php/jmr/article/view/62103

                    %z = ((x-mu)+1i*gamma)/(omega_g);
                    z = ((x-mu)+1i*gamma)/(sigma*sqrt(2));
                    v = (1/(sigma*sqrt(2*pi))) * real(fadf(z)); % Get Voigt from Faddeeva fn.
                case 'num'
                    error('not yet implemented')
            end
        end
    case  1
        %if strcmp(method_type,'fadd')
        %    warning('derivative method uses fadd, method option ignorred')
        %end
        if ~any(strcmp(norm_type,{'amp','int'}))
            error('derivatives are only compatable with amp or int normalization methods')
        end
        % not 100% about this formula, the one listed on wikipedia has  z = (xdiff+1i*gamma)/(omega_g*sqrt(2));
        % and does not have the factor of 2 in front of v
        % but that did not agree with my numerical calculation
        % this agrees to 1e-4
        z = ((x-mu)+1i*gamma)/(sigma*sqrt(2));
        v=- (1./((sigma.^2).*sqrt(pi))) .* real(z.*fadf(z));
        
        
%         xdiff=x-mu;
%         omega_g=omega_g*sqrt(2);
%         xdiff=xdiff*sqrt(2);
%         z = (xdiff+1i*gamma)/(omega_g*sqrt(2));
%         v=- (xdiff./((omega_g.^3).*sqrt(2*pi))) .* real(z.*fadf(z))...
%             +(gamma./((omega_g.^3).*sqrt(2*pi))) .* imag(z.*fadf(z));
    case  2
        if method_type~='fadd'
            warning('derivative method uses fadd, method option ignorred')
        end
        if ~any(strcmp(norm_type,{'amp','int'}))
            error('derivatives are only compatable with amp or int normalization methods')
        end
        %error('not yet implmeneted for 2nd derivative')
        xdiff=x-mu;
        z = ((x-mu)+1i*gamma)/(sigma*sqrt(2));
        v= ((xdiff.^2-gamma^2-sigma^2)./(sigma^4)) * (1/(sigma*sqrt(2*pi))) .* real(fadf(z))...
            - (2*xdiff.*gamma/(sigma^4))*(1/(sigma*sqrt(2*pi))) .* imag(fadf(z))...
            + gamma/((sigma^4)*pi);
        v=v;
end

v=col_vec(v);
switch norm_type
    case 'amp'
        % find the amplitude by evaluating the 'fadd' method at zero
        peak_val=voigt_function_1d(0,sigma,gamma,0,1,0,'norm','int','method','fadd');
        v=v/peak_val;
    case 'sum'
        v=v/sum(v(:));
    case 'int'
         v=v;
end
v=v*amp +offset;

end


function check_input_size(x,var_to_check)
    if numel(var_to_check)~=1 && size(var_to_check)~=size(x)
        error('wrong size input must be the same size as x or scalar')
    end
end