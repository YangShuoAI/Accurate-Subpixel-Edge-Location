function statCircle(xCenter, yCenter, radius, ...
    innerIntensity, outerIntensity, e)

%STATCIRCLE Computes error measure about the edge detection on a circle
    
fprintf ('Error measures in a circle of radius %.2f\n', radius);

error = abs(hypot(e.x-xCenter, e.y-yCenter) - radius);
fprintf ('   Position error:    max=%.4f;  mean=%.4f;  std=%.4f\n', ...
    max(error), mean(error), std(error));
if (outerIntensity > innerIntensity)
    error = 180/pi * abs(atan2(e.ny,e.nx)-atan2(e.y-yCenter,e.x-xCenter));
else
    error = 180/pi * abs(atan2(e.ny,e.nx)-atan2(yCenter-e.y,xCenter-e.x));
end
bigError = error>180;
error(bigError) = 360-error(bigError);
fprintf ('   Orientation error: max=%.4f;  mean=%.4f;  std=%.4f\n', ...
        max(error), mean(error), std(error));
if (outerIntensity > innerIntensity)
    error = abs((e.i1-e.i0)-(outerIntensity-innerIntensity));
else
    error = abs((e.i1-e.i0)-(innerIntensity-outerIntensity));
end
fprintf ('   Intensity error:   max=%.4f;  mean=%.4f;  std=%.4f\n', ...
        max(error), mean(error), std(error));
if (outerIntensity > innerIntensity)
    error = abs(e.curv-1/radius);
else
    error = abs(-e.curv-1/radius);
end
fprintf ('   Curvature error:   max=%.4f;  mean=%.4f;  std=%.4f\n', ...
        max(error), mean(error), std(error));
r = 1./e.curv;
fprintf ('   Radius value:      max=%.3f;  mean=%.3f;  min=%.3f\n', ...
    max(r), mean(r), min(r));  
fprintf ('\n');
end

