function masked=mask_square(in,limits,invert)
%simple masking function that gets used a lot in our data processing
%  limits  =[[tmin,tmax];[xmin,xmax];[ymin,ymax]] (in seconds and meters)
%  invert  - to flip this to an exclusion region instead of a include
%TODO
% could get a speedup buy checking if any of the limits are inf and then not doing that check

if nargin<3
    invert = 0;
end

counts_size=size(in);
if counts_size(2)~=3
    error('wrong txy size')
end 
if ~isequal(size(limits),[3,2])
    error('wrong txylim size')
end

mask=in(:,1)>limits(1,1) & in(:,1)<limits(1,2) &...
     in(:,2)>limits(2,1) & in(:,2)<limits(2,2) &...
     in(:,3)>limits(3,1) & in(:,3)<limits(3,2);
if invert
    mask = ~mask;
end
masked=in(mask,:);
            
end