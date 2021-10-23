%%test_array_to_rgb


test_img=imread('saturn.png');
sfigure(1)
imshow(test_img)
test_img=im2gray(test_img);
test_img=double(test_img);

cmap=viridis;
rgb_out=array_to_rgb(test_img,colormap,[0,200]);
sfigure(2)
imshow(rgb_out)




%%