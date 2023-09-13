function i = ring (xSize, ySize, xCenter, yCenter, ...
    innerRadius, outerRadius, ...
    innerIntensity, outerIntensity, gridResolution)

[x, y] = meshgrid(0.5:xSize+0.5, 0.5:ySize+0.5);
i = -ones(xSize, ySize);

% compute squared distance of each pixel to the image center
dist2 = (x-xCenter).^2 + (y-yCenter).^2;

% compute pixels completely outside 
or2 = outerRadius^2;
c = (dist2(1:xSize, 1:ySize) > or2) + ...
    (dist2(2:end, 1:ySize) > or2) + ...
    (dist2(1:xSize, 2:end) > or2) + ...
    (dist2(2:end, 2:end) > or2);
i(c==4) = outerIntensity;

% compute pixels completely inside 
ir2 = innerRadius^2;
c = (dist2(1:xSize, 1:ySize) < ir2) + ...
    (dist2(2:end, 1:ySize) < ir2) + ...
    (dist2(1:xSize, 2:end) < ir2) + ...
    (dist2(2:end, 2:end) < ir2);
i(c==4) = outerIntensity;

% compute pixels completely in the ring 
c = (dist2(1:xSize, 1:ySize) < or2 & dist2(1:xSize, 1:ySize) > ir2) + ...
    (dist2(2:end, 1:ySize) < or2 & dist2(2:end, 1:ySize) > ir2) + ...
    (dist2(1:xSize, 2:end) < or2 & dist2(1:xSize, 2:end) > ir2) + ...
    (dist2(2:end, 2:end) < or2 & dist2(2:end, 2:end) > ir2);
i(c==4) = innerIntensity;

% compute contour pixels
delta = 1 / (gridResolution-1);
[dx, dy] = meshgrid(-0.5:delta:0.5, -0.5:delta:0.5);
[x, y] = meshgrid(1:xSize, 1:ySize);
i(i<0) = outerIntensity + (innerIntensity-outerIntensity)* ...
    ringGrid(x(i<0), y(i<0), ir2, or2, xCenter, yCenter, dx, dy);

% show image
%imshow(uint8(i));
%axis('image');
end