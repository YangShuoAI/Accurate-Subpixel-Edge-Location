function bw = subpixelImage(edges, dim, factor)
%SUBPIXELIMAGE computes a high resolution edge binary image
%
%   BW = SUBPIXELIMAGE(EDGES, DIM, FACTOR) computes a high resolution edge
%   binary image using the detected EDGES, where DIM is the original image
%   size and FACTOR is the zoom factor
%
%   Example:
%   edges = subpixelEdges(image, 20);
%   bw = subpixelImage(edges, size(image), 5);

rows = factor * dim(1);
columns = factor * dim(2);
bw = logical(zeros(rows, columns));
radii = 1 ./ edges.curv;
radii(radii==inf) = 1e6;
minRadii = 2;
radii(radii>0 & radii<minRadii) = minRadii;
radii(radii<0 & radii>-minRadii) = -minRadii;
center = [edges.x-radii.*edges.nx edges.y-radii.*edges.ny];

for i=1:size(edges.x)
    if abs(edges.ny(i)) > abs(edges.nx(i))
        x = edges.x(i) + linspace(-0.5+0.5/factor,0.5-0.5/factor,factor);
        if center(i,2) > edges.y(i)
            y = center(i,2)-sqrt(radii(i)^2-(x-center(i,1)).^2);
        else
            y = center(i,2)+sqrt(radii(i)^2-(x-center(i,1)).^2);
        end
    else
        y = edges.y(i) + linspace(-0.5+0.5/factor,0.5-0.5/factor,factor);
        if center(i,1) > edges.x(i)
            x = center(i,1)-sqrt(radii(i)^2-(y-center(i,2)).^2);
        else
            x = center(i,1)+sqrt(radii(i)^2-(y-center(i,2)).^2);
        end
    end
    x = (x-0.5)*factor + 0.5;
    y = (y-0.5)*factor + 0.5;
    bw((round(x)-1)*rows+round(y)) = true;
end

