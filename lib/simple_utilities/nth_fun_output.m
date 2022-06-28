function value = nth_fun_output(N,fcn,varargin)
% also see outselect
  [value{1:N}] = fcn(varargin{:});
  value = value{N};
end
