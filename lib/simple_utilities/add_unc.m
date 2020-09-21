function val_unc = add_unc(varargin)
% Accepts tuples of the form [x,x_unc] and returns the pair [val,unc].
% Acts as if all uncertainties are non-negative and adds in quadrature.
%``` Basic function
% x = [0,0.5];
% y = [-1,0.2];
% val_unc = add_unc(x,y);
% isequal(val_unc,[-1, sqrt(0.5^2+0.7^2)])
%```
% Also works for: 
%``` differences
% x = [0,0.5];
% y = [-1,0.2];
% val_unc = add_unc(x,-y);
% isequal(val_unc,[1,sqrt(0.5^2+0.2^2)])
%```
%``` multiple inputs
% x = [0,0.5];
% y = [-1,0.2];
% z = [5,1];
% val_unc = add_unc(x,y,z);
% isequal(val_unc,[4,sqrt(0.5^2+0.2^2+1)])
%```
% Possible improvements: 
% Vectorize to accept two Mx2 arrays or an Nx2 and 1x2 array.
% Implement proper type checking

    if nargin == 1
        warning("Not enough inputs to add_unc");
        if length(varargin{1}) == 2
            val_unc = varargin{1};
        else
            warning("Input to add_unc was not a [v,u] pair");
        end
    elseif nargin == 2
        val_unc = add_unc_pair(varargin{1},varargin{2});
    else
        val_unc = add_unc(varargin{1},add_unc(varargin{2:end}));
    end

    
end

function val_unc = add_unc_pair(x,y)
    val_unc =[x(1) + y(1),sqrt(abs(x(2)).^2+abs(y(2)).^2)];
end
% add_unc = @(x,y) [x(1)+y(1),abs(x(2))+abs(y(2))];