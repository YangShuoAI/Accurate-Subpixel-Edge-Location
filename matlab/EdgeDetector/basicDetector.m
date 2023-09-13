function ep = basicDetector(F, threshold, order)
%BasicDetector computes subpixel edge detection
% F: input image
% threshold: minimum intensity change to be considered as an edge
% order: detect linear edges (1) or second order edges(2)

% initialization
ep = EdgePixel;
[rows, cols] = size(F);
[x, y] = meshgrid(1:cols, 1:rows);

% compute partial derivatives
Fx = zeros(rows, cols);
Fx(1:rows,2:cols-1) = 0.5 *(F(1:rows,3:cols) - F(1:rows,1:cols-2));
Fy = zeros(rows, cols);
Fy(2:rows-1,1:cols) = 0.5 *(F(3:rows,1:cols) - F(1:rows-2,1:cols));
grad = sqrt(Fx.^2+Fy.^2);

% detecte edge pixels with maximum Fy (not including margins)
absFyInner = abs(Fy(3:rows-2,2:cols-1));
E = false(rows, cols);
E(3:rows-2,2:cols-1) = grad(3:rows-2,2:cols-1)>threshold & ...
    absFyInner>=abs(Fx(3:rows-2,2:cols-1)) & ...
    absFyInner>=abs(Fy(2:rows-3,2:cols-1)) & ...
    absFyInner>abs(Fy(4:rows-1,2:cols-1));
edges = (x(E)-1)*rows+y(E);
m = ones(size(edges,1),1);
m(Fx(edges).*Fy(edges)<0) = -1;
SL = zeros(size(edges,1),1);
SM = zeros(size(edges,1),1);
SR = zeros(size(edges,1),1);
for n=-2:2
    SL = SL + F(edges-rows+n);
    SM = SM + F(edges+n);
    SR = SR + F(edges+rows+n);
end
A = (F(edges+2)+F(edges+m*rows+2)+F(edges+m*rows+1)) / 3;
B = (F(edges-m*rows-1)+F(edges-m*rows-2)+F(edges-2)) / 3;
den = 2 * (A-B);
if (order==2)
    c = (SL+SR-2*SM) ./ den;
else
    c = 0;
end
b = (SR-SL) ./ den;
a = (2*SM-5*(A+B)) ./ den - c/12;
n = ones(size(edges,1),1);
n(Fy(edges)<0) = -1;
ep.position = [ep.position; edges];
ep.x = [ep.x; x(edges)];
ep.y = [ep.y; y(edges) - a];
ep.nx = [ep.nx; sign(A-B)./sqrt(1+b.^2).*b];
ep.ny = [ep.ny; sign(A-B)./sqrt(1+b.^2)];
ep.curv = [ep.curv; 2*c.*n./((1+b.^2).^1.5)];
ep.i0 = [ep.i0; min(A,B)];
ep.i1 = [ep.i1; max(A,B)];

% detecte edge pixels with maximum Fx (not including margins)
absFxInner = abs(Fx(2:rows-1,3:cols-2));
E = false(rows, cols);
E(2:rows-1,3:cols-2) = grad(2:rows-1,3:cols-2)>threshold & ...
    absFxInner>abs(Fy(2:rows-1,3:cols-2)) & ...
    absFxInner>=abs(Fx(2:rows-1,2:cols-3)) & ...
    absFxInner>abs(Fx(2:rows-1,4:cols-1));
edges = (x(E)-1)*rows+y(E);
m = ones(size(edges,1),1);
m(Fx(edges).*Fy(edges)<0) = -1;
SL = zeros(size(edges,1),1);
SM = zeros(size(edges,1),1);
SR = zeros(size(edges,1),1);
for n=-2:2
    SL = SL + F(edges-1+n*rows);
    SM = SM + F(edges+n*rows);
    SR = SR + F(edges+1+n*rows);
end
A = (F(edges+2*rows)+F(edges+2*rows+m)+F(edges+rows+m)) / 3;
B = (F(edges-rows-m)+F(edges-2*rows-m)+F(edges-2*rows)) / 3;
den = 2 * (A-B);
if (order==2)
    c = (SL+SR-2*SM) ./ den;
else
    c = 0;
end
b = (SR-SL) ./ den;
a = (2*SM-5*(A+B)) ./ den - c/12;
n = ones(size(edges,1),1);
n(Fx(edges)<0) = -1;
ep.position = [ep.position; edges];
ep.x = [ep.x; x(edges) - a];
ep.y = [ep.y; y(edges)];
ep.nx = [ep.nx; sign(A-B)./sqrt(1+b.^2)];
ep.ny = [ep.ny; sign(A-B)./sqrt(1+b.^2).*b];
ep.curv = [ep.curv; 2*c.*n./((1+b.^2).^1.5)];
ep.i0 = [ep.i0; min(A,B)];
ep.i1 = [ep.i1; max(A,B)];

% show edges
%imshow(uint8(F));
%axis('image');
%hold on;
%seg = 0.6;
%quiver(ep.x-seg/2*ep.ny, ep.y+seg/2*ep.nx, seg*ep.ny, -seg*ep.nx, 0, 'r.');
%quiver(ep.x, ep.y, ep.nx, ep.ny, 0, 'b');
%hold off;
end
