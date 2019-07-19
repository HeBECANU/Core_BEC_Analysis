


out_string=string_value_with_unc(50,1,'s');
if ~strcmp(out_string,'50.0±1.0')
    error('test failed')
end

out_string=string_value_with_unc(50.1234,1,'s');
if ~strcmp(out_string,'50.1±1.0')
    error('test failed')
end

out_string=string_value_with_unc(50.1234,1.8235,'s');
if ~strcmp(out_string,'50±2')
    error('test failed')
end

out_string=string_value_with_unc(50.1234,5.1235,'s');
if ~strcmp(out_string,'50±5')
    error('test failed')
end

out_string=string_value_with_unc(50.1234,0.15,'s');
if ~strcmp(out_string,'50.12±0.15')
    error('test failed')
end
out_string=string_value_with_unc(854738.96541,200,'s');
if ~strcmp(out_string,'854700±200')
    error('test failed')
end

out_string=string_value_with_unc(8.547896541,0.001,'b');
if ~strcmp(out_string,'8.5479(10)')
    error('test failed')
end

out_string=string_value_with_unc(854738.96541,200,'b');
if ~strcmp(out_string,'854700(200)')
    error('test failed')
end

out_string=string_value_with_unc(8.547896541,0.005,'b');
if ~strcmp(out_string,'8.548(5)')
    error('test failed')
end

out_string=string_value_with_unc(854738.96541,130,'b');
if ~strcmp(out_string,'854740(130)')
    error('test failed')
end