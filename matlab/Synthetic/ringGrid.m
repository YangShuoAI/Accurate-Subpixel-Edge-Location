function p = ringGrid(x, y, innerRadius2, outerRadius2, ...
    xCenter, yCenter, dx, dy)
numPixels = size(x, 1);
if (numPixels > 0)
    
    for n=1:numPixels
        value2 = (x(n)+dx-xCenter).^2 + (y(n)+dy-yCenter).^2;
        grid = (value2 < outerRadius2) & (value2 > innerRadius2);
        p(n) = mean(mean(grid));
    end
    
    
    %{
    DX = gpuArray(dx);
    DY = gpuArray(dy);
        p = arrayfun(@(x,y) mean(mean((x+DX-xCenter).^2 + (y+DY-yCenter).^2 < outerRadius2 & ...
        (x+DX-xCenter).^2 + (y+DY-yCenter).^2 > innerRadius2)), x, y);
    %}

    
else
    p = [];
end
end

