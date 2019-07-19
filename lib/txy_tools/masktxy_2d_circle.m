function txymasked=masktxy_2d_circle(txyin,mask_xyr_negpos)
% circular masking function with multiple block allow
% inputs
% txyin - [nx2] vector for time,x,y
% mask_xyr_allow  - [nx4] ech row has xcen,ycen,radius,posneg
% negpos - a true allows everything within radius of xcen,ycen whereas false
% blocks everything within this (but allows all others) 


txy_size=size(txyin);
if txy_size(2)~=3
    error('wrong txy size')
end 
if ~isequal(size(mask_xyr_negpos,2),4)
    error('wrong txylim size')
end

mask_all=true(txy_size(1),1);

xy_in=txyin(:,2:3);

iimax=size(mask_xyr_negpos,1);
for ii=1:iimax
    xcen=mask_xyr_negpos(ii,1);
    ycen=mask_xyr_negpos(ii,2);
    radius=mask_xyr_negpos(ii,3);
    mask_polarity=mask_xyr_negpos(ii,4);
    
    xy_tmp=xy_in; 
    xy_tmp(:,1)=xy_tmp(:,1)-xcen;
    xy_tmp(:,2)=xy_tmp(:,2)-ycen;
    radial_diff=vecnorm(xy_tmp,2,2);
    
    single_mask=radial_diff<radius;
    
    if ~mask_polarity
        single_mask=~single_mask;
    end 
    
    mask_all=mask_all& single_mask;
end

txymasked=txyin(mask_all,:);
            
end