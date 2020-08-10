%test_bib_to_json
 
mpath = fileparts(mfilename('fullpath'));

bibs_dir=fullfile(mpath,'test_bibs')
%bibs_dir='.\figs\to_val_history\papers';

a=bib_to_st(bibs_dir)
for ii=1:numel(a)
    a{ii}
end
    




%function split_ignoring_groups