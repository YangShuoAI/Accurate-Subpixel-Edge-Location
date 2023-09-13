function statRamp(xCenter, yCenter, angle, ...
    innerIntensity, outerIntensity, e)

%STATRAMP Computes error measure about the edge detection on a ramp
    
fprintf ('Error measures in a ramp of angle %.1f\n', angle);

% equation of the ramp: ax+by+c=0
rads = angle/180*pi;
a = -sin(rads);
b = cos(rads);
c = -a*xCenter - b*yCenter;

% position error
error = abs(a*e.x+b*e.y+c) / hypot(a,b);
fprintf ('   Position error:    max=%.4f;  mean=%.4f;  std=%.4f\n', ...
    max(error), mean(error), std(error));


if (outerIntensity > innerIntensity)
    error = 180/pi * abs(atan2(e.ny,e.nx)-atan2(e.y-yCenter,e.x-xCenter));
else
    error = 180/pi * abs(atan2(e.ny,e.nx)-atan2(yCenter-e.y,xCenter-e.x));
end

% orientation error
angle1 = atan2(e.ny, e.nx) / pi * 180;
if (outerIntensity > innerIntensity)
    angle2 = angle - 90;
    if (angle2 < -180)
        angle2 = angle2 + 360;
    end
else
    angle2 = angle + 90;
    if (angle2 > 180)
        angle2 = angle2 - 360;
    end
end
error = abs(angle1 - angle2);
fprintf ('   Orientation error: max=%.4f;  mean=%.4f;  std=%.4f\n', ...
    max(error), mean(error), std(error));

% intensity error
error = abs(abs(e.i1-e.i0)-abs(outerIntensity-innerIntensity));
fprintf ('   Intensity error:   max=%.4f;  mean=%.4f;  std=%.4f\n', ...
        max(error), mean(error), std(error));
fprintf ('\n');

end

