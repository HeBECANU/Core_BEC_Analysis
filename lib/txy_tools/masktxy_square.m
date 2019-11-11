function txymasked=masktxy_square(txyin,txylim)
%simple masking function that gets used a lot in our data processing
%txylim=[[tmin,tmax];[xmin,xmax];[ymin,ymax]] (in seconds and meters)

txymasked=mask_square(txyin,txylim,0);
            
end