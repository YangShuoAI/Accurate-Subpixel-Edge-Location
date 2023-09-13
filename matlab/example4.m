%% SUBPIXEL EDGES - EXAMPLE 4 -----------
% CREATION OF A HIGH RESOLUTION BINARY IMAGE STARTING FROM SUBPIXEL EDGE 
% FEATURES DETECTED ON A IMAGE

addpath(genpath('.'));

%% syntethic ring
imageSize = 35;
xCenter = imageSize/2;
yCenter = imageSize/2;
innerRadius = 8.0;
outerRadius = 10.0;
innerIntensity = 100;
outerIntensity = 200;
gridResolution = 100;
image = ring(imageSize, imageSize, xCenter, yCenter, ...
    innerRadius, outerRadius, innerIntensity, outerIntensity, ...
    gridResolution);

%% subpixel detection
threshold = 25;
edges = subpixelEdges(image, threshold, 'SmoothingIter', 0); 
figure(1);
imshow(image/255, 'InitialMagnification', 'fit');

%% show high resolution binary image
bw = subpixelImage(edges, size(image), 10);
figure(2);
imshow(bw, 'InitialMagnification', 'fit');
