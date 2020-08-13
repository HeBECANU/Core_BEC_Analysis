%test_hotspot_mask

num_pts_in=1e6;
txy_max=40e-3;
sim_counts=2*(rand(num_pts_in,3)-0.5)*txy_max;
sim_data=[];
sim_data.counts_txy={sim_counts};
counts_out=hotspot_mask(sim_data);
pass_frac=counts_out.masked.num_counts/num_pts_in;


% in comparison the pass fraction for a circle inscribed in a quare of size txy_max*2
cicle_pass_frac=@(r) pi *r^2 /((2*txy_max).^2);
circ_rad=txy_max;
rel_pass_frac=pass_frac/cicle_pass_frac(circ_rad);

fprintf('rel hotspot mass pass fraction %.4f  (relative to r=%.2f mm det)\n',rel_pass_frac,circ_rad*1e3)
fprintf('or alternately\n')
circ_rad=38e-3;
rel_pass_frac=pass_frac/cicle_pass_frac(38e-3);
fprintf('rel hotspot mass pass fraction %.4f  (relative to r=%.2f mm det)\n',rel_pass_frac,circ_rad*1e3)