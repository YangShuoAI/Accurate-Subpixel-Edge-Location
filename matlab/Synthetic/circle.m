function i = circle (xSize, ySize, xCenter, yCenter, ...
    radius, innerIntensity, outerIntensity, gridResolution)

% compute pixels completely outside or inside
r2 = radius^2;
[x, y] = meshgrid(1:xSize, 1:ySize);
c = ((x-0.5-xCenter).^2 + (y-0.5-yCenter).^2 < r2);
c = c + ((x-0.5-xCenter).^2 + (y+0.5-yCenter).^2 < r2);
c = c + ((x+0.5-xCenter).^2 + (y-0.5-yCenter).^2 < r2);
c = c + ((x+0.5-xCenter).^2 + (y+0.5-yCenter).^2 < r2);
i = c;
i(c==0) = outerIntensity;
i(c==4) = innerIntensity;

% compute contour pixels
delta = 1 / (gridResolution-1);
[dx, dy] = meshgrid(-0.5:delta:0.5, -0.5:delta:0.5);
i(c>0 & c<4) = outerIntensity + (innerIntensity-outerIntensity)* ...
    circleGrid(x(c>0&c<4),y(c>0&c<4),r2,xCenter,yCenter,dx,dy);

% show image
%imshow(uint8(i));
%axis('image');

% display info
fprintf('Generated circle:\n');
fprintf('   Radius=%.2f;  Center=(%.2f, %.2f)\n', ...
    radius, xCenter, yCenter);
fprintf('   Intensities in/out = (%.2f, %.2f)\n', ...
    innerIntensity, outerIntensity);
fprintf('   Pixel grid resolution = %d\n', gridResolution);
fprintf('\n');
end