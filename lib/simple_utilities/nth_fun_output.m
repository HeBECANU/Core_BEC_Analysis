function value = nth_fun_output(N,fcn,varargin)
  [value{1:N}] = fcn(varargin{:});
  value = value{N};
end
