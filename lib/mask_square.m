function masked=mask_square(in,lim,varargin)
%simple masking function that gets used a lot in our data processing
%txylim=[[tmin,tmax];[xmin,xmax];[ymin,ymax]] (in seconds and meters)
if length(varargin)<1
    invert = 0;
else
    invert = varargin{1};
end

counts_size=size(in);
if counts_size(2)~=3
    error('wrong txy size')
end 
if ~isequal(size(lim),[3,2])
    error('wrong txylim size')
end

mask=in(:,1)>lim(1,1) & in(:,1)<lim(1,2) &...
     in(:,2)>lim(2,1) & in(:,2)<lim(2,2) &...
     in(:,3)>lim(3,1) & in(:,3)<lim(3,2);
if invert
    mask = ~mask;
end
masked=in(mask,:);
            
end