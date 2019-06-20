% A simple test for the conversion between TXY and V_ZXY as commonly used.

% User inputs
n_atoms = 1e3;
error_tolerance = 1e-3;
test_num_exponents = 6:2:24;
num_trials = 1e2;
t_evol = 0.417;
std_v = 1e-2;
std_x = 1e-2;

fprintf(' - Beginning test\n')

h_0 = 0.5*9.81*t_evol^2;
X0 = std_x*randn(n_atoms,3)+h_0; % Initial positions
V0 = std_v*randn(n_atoms,3); % Initial velocities
X1 = nan*X0;
X1(:,[2,3]) = V0(:,[2,3])*t_evol; % free evolution in XY
X1(:,1) = (-V0(:,1) + sqrt(V0(:,1).^2 + 2*9.81*h_0))/9.81; % Analytic expression for arrival time
err = mean(norm(V0-txy_to_vel(X1,0,9.81,h_0))); % Perform the test
logic_string = {'PASS','FAIL'};
pass_val = err < error_tolerance;
fprintf('INFO: Mean error %f, tolerance %f\n',err,error_tolerance)
fprintf('TEST: Mean error below threshold: %s\n',logic_string{pass_val})

fprintf(' - Test complete\n')