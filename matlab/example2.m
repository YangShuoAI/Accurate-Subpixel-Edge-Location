%% SUBPIXEL EDGES - EXAMPLE 2 -----------
% SUBPIXEL EDGE DETECTION IN A REAL IMAGE

addpath(genpath('.'));

%% load image
url='http://up.vbiran.ir/uploads/136213924137813_coins.png';
image = imread(url);
if numel(size(image))==3
    image = rg2gray(image);
end
imshow(image, 'InitialMagnification', 'fit');

%% subpixel detection
threshold = 25;
edges = subpixelEdges(image, threshold); 

%% show edges
visEdges(edges);
