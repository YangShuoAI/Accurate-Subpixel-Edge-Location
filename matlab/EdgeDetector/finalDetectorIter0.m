function ep = finalDetectorIter0(F, threshold, order)
%FinalDetectorIter0 computes subpixel edge detection using floating
% limit windows
% F: input image
% threshold: minimum intensity change to be considered as an edge
% order: detect linear edges (1) or second order edges(2)

%% initialization
ep = EdgePixel;
[rows, cols] = size(F);
[x, y] = meshgrid(1:cols, 1:rows);

% compute partial derivatives
Fx = zeros(rows, cols);
Fx(1:rows,2:cols-1) = 0.5 *(F(1:rows,3:cols) - F(1:rows,1:cols-2));
Fy = zeros(rows, cols);
Fy(2:rows-1,1:cols) = 0.5 *(F(3:rows,1:cols) - F(1:rows-2,1:cols));
grad = sqrt(Fx.^2+Fy.^2);

%% detecte edge pixels with maximum Gy (not including margins)
absFyInner = abs(Fy(6:rows-5,3:cols-2));
E = false(rows, cols);
E(6:rows-5,3:cols-2) = grad(6:rows-5,3:cols-2)>threshold & ...
    absFyInner>=abs(Fx(6:rows-5,3:cols-2)) & ...
    absFyInner>=abs(Fy(5:rows-6,3:cols-2)) & ...
    absFyInner>abs(Fy(7:rows-4,3:cols-2));
edges = (x(E)-1)*rows+y(E);
A = zeros(size(edges,1),1);
B = zeros(size(edges,1),1);
a = zeros(size(edges,1),1);
b = zeros(size(edges,1),1);
c = zeros(size(edges,1),1);

% compute all horizontal edges
for k=1:size(edges,1)
    edge = edges(k);

    % compute window floating limits
    m1=-1;  m2=1;
    if Fx(edge)*Fy(edge)>=0
        m = 1;
        l1=0;       r2=0;
        minl1=-3;   maxr2=3;
        l2=1;       r1=-1;
        maxl2=4;    minr1=-4;
    else
        m = -1;
        l1=-1;      r2=1;
        minl1=-4;   maxr2=4;
        l2=0;       r1=0;
        maxl2=3;    minr1=-3;
    end        
    if abs(Fx(edge))<1
        l1=-1; l2=1; r1=-1; r2=1;
    end
    while l1>minl1 && abs(Fy(edge-rows+l1))>=abs(Fy(edge-rows+l1-1))  
        l1 = l1 - 1;
    end
    while l2<maxl2 && abs(Fy(edge-rows+l2))>=abs(Fy(edge-rows+l2+1))  
        l2 = l2 + 1;
    end
    while m1>-4 && abs(Fy(edge+m1))>=abs(Fy(edge+m1-1))
        m1 = m1 - 1;
    end
    while m2<4 && abs(Fy(edge+m2))>=abs(Fy(edge+m2+1))
        m2 = m2 + 1;
    end
    while r1>minr1 && abs(Fy(edge+rows+r1))>=abs(Fy(edge+rows+r1-1))
        r1 = r1 - 1;
    end
    while r2<maxr2 && abs(Fy(edge+rows+r2))>=abs(Fy(edge+rows+r2+1))
        r2 = r2 + 1;
    end
    
    % compute intensities
    if m > 0
        AA = (F(edge+m2) + F(edge+rows+r2)) / 2;
        BB = (F(edge-rows+l1) + F(edge+m1)) / 2;
    else
        AA = (F(edge-rows+l2) + F(edge+m2)) / 2;
        BB = (F(edge+m1) + F(edge+rows+r1)) / 2;
    end

    % sum columns
    SL=0; SM=0; SR=0;
    for n=l1:l2 SL = SL + F(edge-rows+n); end;
    for n=m1:m2 SM = SM + F(edge+n); end;
    for n=r1:r2 SR = SR + F(edge+rows+n); end;
    
    % compute edge features
    den = 2 * (AA-BB);
    if order == 2 
        c(k) = (SL+SR-2*SM + AA*(2*m2-l2-r2) - BB*(2*m1-l1-r1)) / den;
    else
        c(k) = 0; 
    end
    b(k) = (SR-SL + AA*(l2-r2) - BB*(l1-r1)) / den;
    a(k) = (2*SM - AA*(1+2*m2) - BB*(1-2*m1)) / den - c(k)/12;
    A(k)=AA; B(k)=BB;
    
    
%      if x(edge)==11 && y(edge)==7
%          fprintf ('a=%f b=%f c=%f\n',a(k),b(k),c(k));
%          fprintf('pixel central: Fx=%f  Fy=%f\n', Fx(edge), Fy(edge));
%          fprintf('pixel derecho: Fx=%f  Fy=%f\n', Fx(edge+rows), Fy(edge+rows));
%      end
   
end

% save edges info
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


%% detecte edge pixels with maximum Gx (not including margins)
absFxInner = abs(Fx(3:rows-2,6:cols-5));
E = false(rows, cols);
E(3:rows-2,6:cols-5) = grad(3:rows-2,6:cols-5)>threshold & ...
    absFxInner>abs(Fy(3:rows-2,6:cols-5)) & ...
    absFxInner>=abs(Fx(3:rows-2,5:cols-6)) & ...
    absFxInner>abs(Fx(3:rows-2,7:cols-4));
edges = (x(E)-1)*rows+y(E);
A = zeros(size(edges,1),1);
B = zeros(size(edges,1),1);
a = zeros(size(edges,1),1);
b = zeros(size(edges,1),1);
c = zeros(size(edges,1),1);

% compute all horizontal edges
for k=1:size(edges,1)
    edge = edges(k);
 
    % compute window floating limits
    m1=-1;  m2=1;
    if Fx(edge)*Fy(edge)>=0
        m = 1;
        l1=0;       r2=0;
        minl1=-3;   maxr2=3;
        l2=1;       r1=-1;
        maxl2=4;    minr1=-4;
    else
        m = -1;
        l1=-1;      r2=1;
        minl1=-4;   maxr2=4;
        l2=0;       r1=0;
        maxl2=3;    minr1=-3;
    end       
    if abs(Fy(edge))<1
        l1=-1; l2=1; r1=-1; r2=1;
    end
    while l1>minl1 && abs(Fx(edge-1+l1*rows))>=abs(Fx(edge-1+(l1-1)*rows))
        l1 = l1 - 1;
    end
    while l2<maxl2 && abs(Fx(edge-1+l2*rows))>=abs(Fx(edge-1+(l2+1)*rows))  
        l2 = l2 + 1;
    end
    while m1>-4 && abs(Fx(edge+m1*rows))>=abs(Fx(edge+(m1-1)*rows))
        m1 = m1 - 1;
    end
    while m2<4 && abs(Fx(edge+m2*rows))>=abs(Fx(edge+(m2+1)*rows))
        m2 = m2 + 1;
    end
    while r1>minr1 && abs(Fx(edge+1+r1*rows))>=abs(Fx(edge+1+(r1-1)*rows))
        r1 = r1 - 1;
    end
    while r2<maxr2 && abs(Fx(edge+1+r2*rows))>=abs(Fx(edge+1+(r2+1)*rows))
        r2 = r2 + 1;
    end

    % compute intensities
    if m > 0
        AA = (F(edge+m2*rows) + F(edge+1+r2*rows)) / 2;
        BB = (F(edge-1+l1*rows) + F(edge+m1*rows)) / 2;
    else
        AA = (F(edge-1+l2*rows) + F(edge+m2*rows)) / 2;
        BB = (F(edge+m1*rows) + F(edge+1+r1*rows)) / 2;
    end

    % sum rows
    SL=0; SM=0; SR=0;
    for n=l1:l2 SL = SL + F(edge-1+n*rows); end;
    for n=m1:m2 SM = SM + F(edge+n*rows); end;
    for n=r1:r2 SR = SR + F(edge+1+n*rows); end;
   
    % compute edge features
    den = 2 * (AA-BB);
    if order == 2
        c(k) = (SL+SR-2*SM + AA*(2*m2-l2-r2) - BB*(2*m1-l1-r1)) / den;
    else
        c(k) = 0; 
    end
    b(k) = (SR-SL + AA*(l2-r2) - BB*(l1-r1)) / den;
    a(k) = (2*SM - AA*(1+2*m2) - BB*(1-2*m1)) / den - c(k)/12;
    A(k)=AA; B(k)=BB;
end

% save edges info
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
end