function ep = smoothDetector(F, threshold, order)
%SMOOTHDETECTOR computes subpixel edge detection
% F: input image
% threshold: minimum intensity change to be considered as an edge
% degree: detect linear edges (1) or second order edges(2)

% initialization
ep = EdgePixel;
[rows, cols] = size(F);
[x, y] = meshgrid(1:cols, 1:rows);

% smooth image
G = double(F);
G(2:rows-1,2:cols-1) = ...
    (F(1:rows-2,1:cols-2) + F(1:rows-2,2:cols-1) + F(1:rows-2,3:cols) + ...
    F(2:rows-1,1:cols-2) + F(2:rows-1,2:cols-1) + F(2:rows-1,3:cols) + ...
    F(3:rows,1:cols-2) + F(3:rows,2:cols-1) + F(3:rows,3:cols))/9;

% compute partial derivatives
Gx = zeros(rows, cols);
Gx(1:rows,2:cols-1) = 0.5 *(G(1:rows,3:cols) - G(1:rows,1:cols-2));
Gy = zeros(rows, cols);
Gy(2:rows-1,1:cols) = 0.5 *(G(3:rows,1:cols) - G(1:rows-2,1:cols));
grad = sqrt(Gx.^2+Gy.^2);

% detecte edge pixels with maximum Gy (not including margins)
absGyInner = abs(Gy(5:rows-4,2:cols-1));
E = false(rows, cols);
E(5:rows-4,2:cols-1) = grad(5:rows-4,2:cols-1)>threshold & ...
    absGyInner>=abs(Gx(5:rows-4,2:cols-1)) & ...
    absGyInner>abs(Gy(4:rows-5,2:cols-1)) & ...
    absGyInner>abs(Gy(6:rows-3,2:cols-1));
edges = (x(E)-1)*rows+y(E);
m = ones(size(edges,1),1);
m(Gx(edges).*Gy(edges)<0) = -1;
SL = zeros(size(edges,1),1);
SM = zeros(size(edges,1),1);
SR = zeros(size(edges,1),1);
for n=-3:3
    SL = SL + G(edges-rows+n+m);
    SM = SM + G(edges+n);
    SR = SR + G(edges+rows+n-m);
end
A = (G(edges+4)+G(edges+m*rows+4)+G(edges+m*rows+3)) / 3;
B = (G(edges-m*rows-3)+G(edges-m*rows-4)+G(edges-4)) / 3;
den = 2 * (A-B);
if (order==2)
    c = (SL+SR-2*SM) ./ den;
else
    c = 0;
end
b = m + (SR-SL) ./ den;
a = (2*SM-7*(A+B)) ./ den - 3*c/4;    
n = ones(size(edges,1),1);
n(Gy(edges)<0) = -1;
ep.position = [ep.position; edges];
ep.x = [ep.x; x(edges)];
ep.y = [ep.y; y(edges) - a];
ep.nx = [ep.nx; sign(A-B)./sqrt(1+b.^2).*b];
ep.ny = [ep.ny; sign(A-B)./sqrt(1+b.^2)];
ep.curv = [ep.curv; 2*c.*n./((1+b.^2).^1.5)];
ep.i0 = [ep.i0; min(A,B)];
ep.i1 = [ep.i1; max(A,B)];

%{
for k=1:size(edges,1)
    edge = edges(k);    
    if x(edge)==172 && y(edge)==143
        fprintf ('a=%f b=%f c=%f\n',a(k),b(k),c);
    end
end
%}

% detecte edge pixels with maximum Gx (not including margins)
absGxInner = abs(Gx(2:rows-1,5:cols-4));
E = false(rows, cols);
E(2:rows-1,5:cols-4) = grad(2:rows-1,5:cols-4)>threshold & ...
    absGxInner>=abs(Gy(2:rows-1,5:cols-4)) & ...
    absGxInner>abs(Gx(2:rows-1,4:cols-5)) & ...
    absGxInner>abs(Gx(2:rows-1,6:cols-3));
edges = (x(E)-1)*rows+y(E);
m = ones(size(edges,1),1);
m(Gx(edges).*Gy(edges)<0) = -1;
SL = zeros(size(edges,1),1);
SM = zeros(size(edges,1),1);
SR = zeros(size(edges,1),1);
for n=-3:3
    SL = SL + G(edges-1+(n+m)*rows);
    SM = SM + G(edges+n*rows);
    SR = SR + G(edges+1+(n-m)*rows);
end
A = (G(edges+4*rows)+G(edges+4*rows+m)+G(edges+3*rows+m)) / 3;
B = (G(edges-3*rows-m)+G(edges-4*rows-m)+G(edges-4*rows)) / 3;
den = 2 * (A-B);
if (order==2)
    c = (SL+SR-2*SM) ./ den;
else
    c = 0;
end
b = m + (SR-SL) ./ den;
a = (2*SM-7*(A+B)) ./ den - 3*c/4;
n = ones(size(edges,1),1);
n(Gx(edges)<0) = -1;
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

