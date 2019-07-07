function out_string=string_value_with_unc(value,unc,type,remove_common_factors)
%string_value_with_unc - prints value with uncert. and the right number of sig figs
%
% Syntax:  output_string=import_data(value,unc,type)
%
% Inputs:
%    value      - numeric, value 
%    unc        - nmeric, uncertianty in the value
%    type       - to use the ± notation as in 10±3 pass 's','standard','pm' to use the () notation as in 10(3) pass , 'b','brackets'

%
% Outputs:
%    out_string - char vec, the formated value

%
% Example: 

% string_value_with_unc(50.1234,0.15,'s')
% '50.12±0.15'
% string_value_with_unc(8.547896541,0.001,'b')
% '8.5479(10)'
% string_value_with_unc(854738.96541,200,'b')
% '854700(200)'

% Other m-files required: none
% Also See:test_string_value_with_unc
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%    -more commenting
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2019-07-02

%------------- BEGIN CODE --------------


if nargin<3
    type='s';
end

if ~isnumeric(value) || ~isnumeric(unc)
    error('value and uncert. must be numeric')
end
% find the first decimal place
unc_first_decimal_place=floor(log10(unc));
% i think this section can be done a little better
fractional_unc_error=abs(((10^unc_first_decimal_place)-unc)/unc);
if fractional_unc_error<=(1/3)
    two_sig_figs=1;
else
    two_sig_figs=0;
end
% if the value is beteen 1*10^i and 1.5*10^i (where i is an integer) then show two decimal places
decimal_place_unc=unc_first_decimal_place-two_sig_figs;
rounded_unc=(10^decimal_place_unc)*round(unc*(10^(-decimal_place_unc)));
rounded_val=(10^decimal_place_unc)*round(value*(10^(-decimal_place_unc)));

% standard plus-minus format
if sum(strcmp(type,{'s','standard','pm'}))>0

    unc_str=sprintf(cat(2,'%.',sprintf('%u',max([0,-decimal_place_unc])),'f'),rounded_unc);
    val_str=sprintf(cat(2,'%.',sprintf('%u',max([0,-decimal_place_unc])),'f'),rounded_val);
    out_string=cat(2,val_str,'±',unc_str);

    
%metrology brackets notation    
elseif  sum(strcmp(type,{'b','brackets'}))>0
    unc_str=sprintf('%.0f',rounded_unc*(10^(max([0,-decimal_place_unc]))));
    val_str=sprintf(cat(2,'%.',sprintf('%u',max([0,-decimal_place_unc])),'f'),rounded_val);
    out_string=cat(2,val_str,'(',unc_str,')');
else
    error('did not pass valid format string')
end


end