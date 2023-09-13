function fullStatCircle (xSize, ySize, xCenter, yCenter, ...
    radius, innerIntensity, outerIntensity, gridResolution, ...
    noiseMag, threshold, method, order)
%STATCIRCLE Computes error measure about the edge detection on a circle

% create image
i = circle(xSize, ySize, xCenter, yCenter, ...
    radius, innerIntensity, outerIntensity, gridResolution);

% add noise
if (noiseMag > 0)
    i = noise(i, noiseMag);
end

% detect edges
switch method
    case 1
        e = basicDetector(i, threshold, order);
    case 2
        e = smoothDetector(i, threshold, order);
end
    
% error measure
statCircle(xCenter, yCenter, radius, innerIntensity, outerIntensity, e);

end

