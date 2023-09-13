function ep = nearbyEdgeDetector(F, threshold, order)
%NEARBYEDGEDETECTOR computes subpixel edge detection using floating windows
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
    absGyInner>=abs(Gy(4:rows-5,2:cols-1)) & ...
    absGyInner>abs(Gy(6:rows-3,2:cols-1));
edges = (x(E)-1)*rows+y(E);
m = ones(size(edges,1),1);
m(Gx(edges).*Gy(edges)<0) = -1;
minl1 = -3 * ones(size(edges,1),1);
minl1(m<0) = -4;
maxr2 = -minl1;
minr1 = -4 * ones(size(edges,1),1);
minr1(m<0) = -3;
maxl2 = -minr1;
minl2 = 2 * ones(size(edges,1),1);
minl2(m<0) = 1;
maxr1 = -minl2;
maxl1 = -1 * ones(size(edges,1),1);
maxl1(m<0) = -2;
minr2 = -maxl1;
A = zeros(size(edges,1),1);
B = zeros(size(edges,1),1);
a = zeros(size(edges,1),1);
b = zeros(size(edges,1),1);
c = zeros(size(edges,1),1);

% compute all horizontal edges
for k=1:size(edges,1)
    edge = edges(k);
    l1 = maxl1(k);
    while l1 > minl1(k)   
        l1 = l1 - 1;
        if abs(Gy(edge-rows+l1))>=abs(Gy(edge-rows)) || ...
                sign(Gy(edge-rows+l1)*Gy(edge-rows))==-1
            l1 = l1 + 1;
            break;
        end
    end
    l2 = minl2(k);
    while l2 < maxl2(k)   
        l2 = l2 + 1;
        if abs(Gy(edge-rows+l2))>=abs(Gy(edge-rows)) || ...
                sign(Gy(edge-rows+l2)*Gy(edge-rows))==-1
            l2 = l2 - 1;
            break;
        end
    end
    m1 = -1;
    while m1 > -4   
        m1 = m1 - 1;
        if abs(Gy(edge+m1))>=abs(Gy(edge)) || ...
                sign(Gy(edge+m1)*Gy(edge))==-1
            m1 = m1 + 1;
            break;
        end
    end
    m2 = 1;
    while m2 < 4   
        m2 = m2 + 1;
        if abs(Gy(edge+m2))>=abs(Gy(edge)) || ...
                sign(Gy(edge+m2)*Gy(edge))==-1
            m2 = m2 - 1;
            break;
        end
    end
    r1 = maxr1(k);
    while r1 > minr1(k)   
        r1 = r1 - 1;
        if abs(Gy(edge+rows+r1))>=abs(Gy(edge+rows)) || ...
                sign(Gy(edge+rows+r1)*Gy(edge+rows))==-1
            r1 = r1 + 1;
            break;
        end
    end
    r2 = minr2(k);
    while r2 < maxr2(k)   
        r2 = r2 + 1;
        if abs(Gy(edge+rows+r2))>=abs(Gy(edge+rows)) || ...
                sign(Gy(edge+rows+r2)*Gy(edge+rows))==-1
            r2 = r2 - 1;
            break;
        end
    end
    SL=0; SM=0; SR=0;
    for n=l1:l2 SL = SL + G(edge-rows+n); end;
    for n=m1:m2 SM = SM + G(edge+n); end;
    for n=r1:r2 SR = SR + G(edge+rows+n); end;
    if m(k) > 0
        AA = (G(edge+m2) + G(edge+rows+r2)) / 2;
        BB = (G(edge-rows+l1) + G(edge+m1)) / 2;
    else
        AA = (G(edge-rows+l2) + G(edge+m2)) / 2;
        BB = (G(edge+m1) + G(edge+rows+r1)) / 2;
    end
    den = 2 * (AA-BB);
    if (order==2) 
        c(k) = (SL+SR-2*SM + AA*(2*m2-l2-r2) - BB*(2*m1-l1-r1)) / den;
    else
        c(k) = 0; 
    end
    b(k) = (SR-SL + AA*(l2-r2) - BB*(l1-r1)) / den;
    a(k) = (2*SM - AA*(1+2*m2) - BB*(1-2*m1)) / den - 3*c(k)/4;
    A(k)=AA; B(k)=BB;
    
    %{
    if x(edge)==172 && y(edge)==143
        fprintf ('a=%f b=%f c=%f\n',a,b,c);
    end
   %}
end

% save edges info
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


% detecte edge pixels with maximum Gx (not including margins)
absGxInner = abs(Gx(2:rows-1,5:cols-4));
E = false(rows, cols);
E(2:rows-1,5:cols-4) = grad(2:rows-1,5:cols-4)>threshold & ...
    absGxInner>abs(Gy(2:rows-1,5:cols-4)) & ...
    absGxInner>=abs(Gx(2:rows-1,4:cols-5)) & ...
    absGxInner>abs(Gx(2:rows-1,6:cols-3));
edges = (x(E)-1)*rows+y(E);
m = ones(size(edges,1),1);
m(Gx(edges).*Gy(edges)<0) = -1;
minl1 = -3 * ones(size(edges,1),1);
minl1(m<0) = -4;
maxr2 = -minl1;
minr1 = -4 * ones(size(edges,1),1);
minr1(m<0) = -3;
maxl2 = -minr1;
minl2 = 2 * ones(size(edges,1),1);
minl2(m<0) = 1;
maxr1 = -minl2;
maxl1 = -1 * ones(size(edges,1),1);
maxl1(m<0) = -2;
minr2 = -maxl1;
A = zeros(size(edges,1),1);
B = zeros(size(edges,1),1);
a = zeros(size(edges,1),1);
b = zeros(size(edges,1),1);
c = zeros(size(edges,1),1);

% compute all horizontal edges
for k=1:size(edges,1)
    edge = edges(k);
    l1 = maxl1(k);
    while l1 > minl1(k)   
        l1 = l1 - 1;
        if abs(Gx(edge-1+l1*rows))>=abs(Gx(edge-1)) || ...
                sign(Gx(edge-1+l1*rows)*Gx(edge-1))==-1
            l1 = l1 + 1;
            break;
        end
    end
    l2 = minl2(k);
    while l2 < maxl2(k)   
        l2 = l2 + 1;
        if abs(Gx(edge-1+l2*rows))>=abs(Gx(edge-1)) || ...
                sign(Gx(edge-1+l2*rows)*Gx(edge-1))==-1
            l2 = l2 - 1;
            break;
        end
    end
    m1 = -1;
    while m1 > -4   
        m1 = m1 - 1;
        if abs(Gx(edge+m1*rows))>=abs(Gx(edge)) || ...
                sign(Gx(edge+m1*rows)*Gx(edge))==-1
            m1 = m1 + 1;
            break;
        end
    end
    m2 = 1;
    while m2 < 4   
        m2 = m2 + 1;
        if abs(Gx(edge+m2*rows))>=abs(Gx(edge)) || ...
                sign(Gx(edge+m2*rows)*Gx(edge))==-1
            m2 = m2 - 1;
            break;
        end
    end
    r1 = maxr1(k);
    while r1 > minr1(k)   
        r1 = r1 - 1;
        if abs(Gx(edge+1+r1*rows))>=abs(Gx(edge+1)) || ...
                sign(Gx(edge+1+r1*rows)*Gx(edge+1))==-1
            r1 = r1 + 1;
            break;
        end
    end
    r2 = minr2(k);
    while r2 < maxr2(k)   
        r2 = r2 + 1;
        if abs(Gx(edge+1+r2*rows))>=abs(Gx(edge+1)) || ...
                sign(Gx(edge+1+r2*rows)*Gx(edge+1))==-1
            r2 = r2 - 1;
            break;
        end
    end
    SL=0; SM=0; SR=0;
    for n=l1:l2 SL = SL + G(edge-1+n*rows); end;
    for n=m1:m2 SM = SM + G(edge+n*rows); end;
    for n=r1:r2 SR = SR + G(edge+1+n*rows); end;
    if m(k) > 0
        AA = (G(edge+m2*rows) + G(edge+1+r2*rows)) / 2;
        BB = (G(edge-1+l1*rows) + G(edge+m1*rows)) / 2;
    else
        AA = (G(edge-1+l2*rows) + G(edge+m2*rows)) / 2;
        BB = (G(edge+m1*rows) + G(edge+1+r1*rows)) / 2;
    end
    den = 2 * (AA-BB);
    if (order==2) 
        c(k) = (SL+SR-2*SM + AA*(2*m2-l2-r2) - BB*(2*m1-l1-r1)) / den;
    else
        c(k) = 0; 
    end
    b(k) = (SR-SL + AA*(l2-r2) - BB*(l1-r1)) / den;
    a(k) = (2*SM - AA*(1+2*m2) - BB*(1-2*m1)) / den - 3*c(k)/4;
    A(k)=AA; B(k)=BB;
end

% save edges info
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
end

