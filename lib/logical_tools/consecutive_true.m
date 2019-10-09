function out=consecutive_true(in,n,assum_start_high,assum_end_high)
% given some logical vector 
% return true when the output will go high for n positions
% starting at the first element of that block

if nargin<3
    assum_start_high=false;
end

if nargin<4
    assum_end_high=false;
end
    

in=col_vec(in)'; % to be strfind compatable row vec
% padd the input with zeros on either side
if assum_start_high
    in=[0,true(1,n),in];
else
    in=[0,false(1,n),in];
end

if assum_end_high
    in=[in,true(1,n),0];
else
    in=[in,false(1,n),0];
end

%TODO
%start & end behaviour

%find positions where the logic should transition to true
logic_up=strfind(in,[0,true(1,n)])+1; %+1 to account for the false

logic_down=strfind(in,[true(1,n),0])+n-1;

out=false*in;
if numel(logic_up)~=numel(logic_down)
    error('logic up and down length do not match')
end
for ii=1:numel(logic_up)
    out(logic_up(ii):logic_down(ii))=true;
end
%trim off the inital padding
out=out(n+2:end-n-1);
end