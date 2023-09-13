%% SUBPIXEL EDGES - EXAMPLE 6 -----------
% SUBPIXEL EDGE DETECTION IN A REAL ANGIOGRAPHY

addpath(genpath('.'));

%% load image
% url='http://serdis.dis.ulpgc.es/~atrujillo/ngImgCrop-master/test/angio2.PNG';

url = "angio2.png"; %"image1.png";%"image2.jpeg";%
img_color = imread(url);
image = rgb2gray(img_color);

%% subpixel detection
threshold = 4;
iter = 3;
[edges, RI] = subpixelEdges(image, threshold, 'SmoothingIter', iter); 

%% show image
showRestoredImage = false;
if showRestoredImage
    imshow(RI/255,'InitialMagnification', 'fit');
else
    imshow(img_color,'InitialMagnification', 'fit');
end

%% show edges
visEdges(edges);

