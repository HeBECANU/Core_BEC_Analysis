function lambda=tf_expand_scaling_trap_off_asm_sph_approx(omega_tzero,t)
% long time aproximation for a spherical condensate
% https://iopscience.iop.org/article/10.1088/1674-1056/ab4177
% eq 7
omega_rad=sqrt(mean(omega_tzero.^2));


lambda=sqrt(2/3).*omega_rad.*t + ...
        sqrt(pi) *gamma(2/3) /gamma(1/6);

lambda=repmat(lambda,[1,3]);


end