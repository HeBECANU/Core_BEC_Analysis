


out_string=string_value_with_unc(50,1,'s');
strcmp(out_string,'50±1')




out_string=string_value_with_unc(50.1234,1,'s');
strcmp(out_string,'50±1')




out_string=string_value_with_unc(50.1234,1.8235,'s');
strcmp(out_string,'50±2')


out_string=string_value_with_unc(50.1234,5.1235,'s');
strcmp(out_string,'50±5')

out_string=string_value_with_unc(50.1234,0.15,'s');
strcmp(out_string,'50.12±0.15')

out_string=string_value_with_unc(854738.96541,200,'s')
strcmp(out_string,'854700±200')

out_string=string_value_with_unc(8.547896541,0.001,'b');
strcmp(out_string,'8.5479(10)')

%%
out_string=string_value_with_unc(854738.96541,200,'b')
strcmp(out_string,'854700(200)')


%%

out_string=string_value_with_unc(8.547896541,0.005,'b');
strcmp(out_string,'8.548(5)')


%%

out_string=string_value_with_unc(854738.96541,130,'b')
strcmp(out_string,'854700(200)')
