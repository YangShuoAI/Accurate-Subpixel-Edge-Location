function statRing(xCenter, yCenter, innerRadius, outerRadius, ...
    innerIntensity, outerIntensity, e)

%STATRING Computes error measure about the edge detection on a ring
    
fprintf ('Error measures in a ring of radius %.2f & %.2f\n', ...
    innerRadius, outerRadius);

% position error
dist = hypot(e.x-xCenter, e.y-yCenter);
meanRadius = (innerRadius + outerRadius) / 2;
radius = ones(size(dist)) * innerRadius;
radius(dist > meanRadius) = outerRadius;
error = abs(dist - radius);
fprintf ('   Position error:    max=%.4f;  mean=%.4f;  std=%.4f\n', ...
    max(error), mean(error), std(error));

% orientation error
angle1 = atan2(e.ny, e.nx) / pi * 180;
angle2 = atan2(e.y-yCenter, e.x-xCenter) / pi * 180;
angle2(dist < meanRadius) = angle2(dist < meanRadius) + 180;
if (outerIntensity < innerIntensity)
    angle2 = angle2 + 180;
end
angle2 = mod(angle2+180, 360) - 180;
angle1(angle1==-180) = 180;
angle2(angle2==-180) = 180;
error = abs(angle1 - angle2);
error(error>180) = 360-error(error>180);

fprintf ('   Orientation error: max=%.4f;  mean=%.4f;  std=%.4f\n', ...
    max(error), mean(error), std(error));

% curvature error
if (outerIntensity > innerIntensity)
    error = abs(e.curv-(radius.^-1));
else
    error = abs(-e.curv-(radius.^-1));
end
fprintf ('   Curvature error:   max=%.4f;  mean=%.4f;  std=%.4f\n', ...
        max(error), mean(error), std(error));

% intensity error
error = abs(abs(e.i1-e.i0)-abs(outerIntensity-innerIntensity));
fprintf ('   Intensity error:   max=%.4f;  mean=%.4f;  std=%.4f\n', ...
        max(error), mean(error), std(error));
fprintf ('\n');

end

