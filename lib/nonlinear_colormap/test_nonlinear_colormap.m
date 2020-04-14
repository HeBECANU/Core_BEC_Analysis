% test nonlinear colormap

im_in=imread('spine.tif');


%%
figure(1)
subplot(1,2,1)
imagesc(im_in)
cm=inferno;
colormap(cm)

%%
figure(1)
subplot(1,2,2)
imagesc(im_in)
cm=inferno;
cm_out=nonlinear_colormap(cm,'logistic',[20,0.7,0.5,0,1],1)
figure(1)
subplot(1,2,2)
colormap(gca,cm_out)

%%
figure(1)
subplot(1,2,2)
imagesc(im_in)
cm=inferno;
cm_out=nonlinear_colormap(cm,'logistic',[25,0.8,0.5,nan,1],1)
figure(1)
subplot(1,2,2)
colormap(gca,cm_out)


%%
figure(1)
subplot(1,2,2)
imagesc(im_in)
cm=inferno;
cm_out=nonlinear_colormap(cm,'power',[0.5],1)
figure(1)
subplot(1,2,2)
colormap(gca,cm_out)


