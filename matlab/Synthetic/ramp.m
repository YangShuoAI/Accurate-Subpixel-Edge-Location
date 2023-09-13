function i = ramp(xSize, ySize, xCenter, yCenter, ...
    angle, innerIntensity, outerIntensity)

% equation of the ramp: ax+by+c=0
rads = angle/180*pi;
a = -sin(rads);
b = cos(rads);
c = -a*xCenter - b*yCenter;

% creating image
[x, y] = meshgrid(1:xSize, 1:ySize);
i = a*(x-0.5) + b*(y-0.5) + c > 0;
i = i + (a*(x-0.5) + b*(y+0.5) + c > 0);
i = i + (a*(x+0.5) + b*(y-0.5) + c > 0);
i = i + (a*(x+0.5) + b*(y+0.5) + c > 0);
i(i==0) = outerIntensity;
i(i==4) = innerIntensity;
i(i>0 & i<4) = outerIntensity + (innerIntensity-outerIntensity)* ...
    rampGrid(x(i>0&i<4),y(i>0&i<4),a,b,c,100);

% show image
%imshow(uint8(i));
%axis('image');
end