


string_value_with_unc(50,1)

out_string=string_value_with_unc(50,1,'type','s');
if ~strcmp(out_string,'50.0±1.0')
    error('test failed')
end

out_string=string_value_with_unc(50.1234,1,'type','s');
if ~strcmp(out_string,'50.1±1.0')
    error('test failed')
end

out_string=string_value_with_unc(50.1234,1.8235,'type','s');
if ~strcmp(out_string,'50±2')
    error('test failed')
end

out_string=string_value_with_unc(50.1234,5.1235,'type','s');
if ~strcmp(out_string,'50±5')
    error('test failed')
end

out_string=string_value_with_unc(50.1234,0.15,'type','s');
if ~strcmp(out_string,'50.12±0.15')
    error('test failed')
end
out_string=string_value_with_unc(854738.96541,200,'type','s','separator',',');
if ~strcmp(out_string, '854,700±200')
    error('test failed')
end

out_string=string_value_with_unc(8.547896541,0.001,'type','s');
if ~strcmp(out_string,'8.5479±0.0010')
    error('test failed')
end

out_string=string_value_with_unc(8.547896541,0.001,'type','s','separator',',');
if ~strcmp(out_string,'8.547,9±0.001,0')
    error('test failed')
end

out_string=string_value_with_unc(8.547896541,0.001,'type','b');
if ~strcmp(out_string,'8.547,9(10)')
    error('test failed')
end

out_string=string_value_with_unc(854738.96541,200,'type','b');
if ~strcmp(out_string,'854,700(200)')
    error('test failed')
end

out_string=string_value_with_unc(8.547896541,0.0000005,'type','b');
if ~strcmp(out_string,'8.547,896,5(5)')
    error('test failed')
end

out_string=string_value_with_unc(854738.96541,130,'type','b');
if ~strcmp(out_string,'854,740(130)')
    error('test failed')
end
%% Common factors
fac=1e6;
out_string=string_value_with_unc(854738.96541*fac,130*fac,'type','s','separator',false)

fac=1e6;
out_string=string_value_with_unc(854738.96541*fac,130*fac,'type','b','separator',false)

fac=1e6;
out_string=string_value_with_unc(854738.96541*fac,110*fac,'type','s','separator',false)

