addpath(genpath('../Subpixel Matlab denoising'));

%% syntethic circle
radius = 15.0;
innerIntensity = 100;
outerIntensity = 200;
gridResolution = 100;
imageSize = 2*radius+20;
xCenter = imageSize/2 + 0.15;
yCenter = imageSize/2 + 0.23;
image = circle(imageSize+5, imageSize+5, xCenter, yCenter, ...
    radius, innerIntensity, outerIntensity, gridResolution);
imshow(image/255,'InitialMagnification', 'fit');

threshold = 25;
edges = subpixelEdges(image, threshold, 'SmoothingIter', 0, 'Order', 2); 

visEdges(edges);

statCircle(xCenter, yCenter, radius, ...
    innerIntensity, outerIntensity, edges)