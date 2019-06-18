function txymasked=masktxy(txyin,txylim)
%simple masking function that gets used a lot in our data processing
%txylim=[[tmin,tmax];[xmin,xmax];[ymin,ymax]] (in seconds and meters)
txy_size=size(txyin);
if txy_size(2)~=3
    error('wrong txy size')
end 
if ~isequal(size(txylim),[3,2])
    error('wrong txylim size')
end

mask=txyin(:,1)>txylim(1,1) & txyin(:,1)<txylim(1,2) &...
     txyin(:,2)>txylim(2,1) & txyin(:,2)<txylim(2,2) &...
     txyin(:,3)>txylim(3,1) & txyin(:,3)<txylim(3,2);
txymasked=txyin(mask,:);
            
end