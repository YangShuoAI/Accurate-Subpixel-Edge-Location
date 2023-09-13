function [ep, I] = finalDetectorIterN(F, threshold, order)
%finalDetectorIter applies the iterative restoration process
% It computes subpixel edge detection using floating 
% limit windows, and detecting cases with very close edges.
% In every iteration, a synthetic image is created to be used in the next
% iteration
% F: input image
% threshold: minimum intensity change to be considered as an edge
% degree: detect linear edges (1) or second order edges(2)

% DECIR QUE SE PUEDE LLAMAR COMO [EDGES, IMAGE] = ....


%% initialization
ep = EdgePixel;
[rows, cols] = size(F);
[x, y] = meshgrid(1:cols, 1:rows);
I = zeros(rows, cols);
C = zeros(rows, cols);
pixelGridResol = 50;

% smooth image
w = 0.75;  % case for smooth mask 'average', H(1:3,1:3) = 1/9
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


%% detecte edge pixels with maximum Gy (not including margins)
absGyInner = abs(Gy(6:rows-5,3:cols-2));
E = false(rows, cols);
E(6:rows-5,3:cols-2) = grad(6:rows-5,3:cols-2)>threshold & ...
    absGyInner>=abs(Gx(6:rows-5,3:cols-2)+1e-8) & ...
    absGyInner>=abs(Gy(5:rows-6,3:cols-2)+1e-8) & ...
    absGyInner>abs(Gy(7:rows-4,3:cols-2)+1e-8);

% % E(6:rows-5,3:cols-2) = absGyInner>=abs(Gx(6:rows-5,3:cols-2));

% [rows1, cols1] = size(E);
% fileID = fopen('E.txt','w+');
% for i = 1:rows1
%     for j = 1:cols1
%         fprintf(fileID, '%d ', E(i,j));
%     end
%     fprintf(fileID, '\n');
% end
% fclose(fileID);

% fileID = fopen('grad.txt','w+');
% for i = 1:rows
%     for j = 1:cols
%         fprintf(fileID, '%.4f ', grad(i,j));
%     end
%     fprintf(fileID, '\n');
% end
% fclose(fileID);

flag = x(E);
edges = (x(E)-1)*rows+y(E);
A = zeros(size(edges,1),1);
B = zeros(size(edges,1),1);
a = zeros(size(edges,1),1);
b = zeros(size(edges,1),1);
c = zeros(size(edges,1),1);
nx = zeros(size(edges,1),1);
ny = zeros(size(edges,1),1);
curv = zeros(size(edges,1),1);
valid = logical(zeros(size(edges,1),1));


% fileID = fopen('edges.txt','w+');
% for iii = 1:size(edges,1)
%     fprintf(fileID, '%d\n', edges(iii));
% end
% fclose(fileID);

% compute all horizontal edges
for k=1:size(edges,1)
    if k == 100
        1 == 1;
    end
    edge = edges(k);

    % compute window floating limits
    m1=-1;  m2=1;
    if Gx(edge)*Gy(edge)>=0
        m = 1;
        l1=-1;      r2=1;
        minl1=-3;   maxr2=3;
        l2=1;       r1=-1;
        maxl2=4;    minr1=-4;
    else
        m = -1;
        l1=-1;      r2=1;
        minl1=-4;   maxr2=4;
        l2=1;       r1=-1;
        maxl2=3;    minr1=-3;
    end
    while l1>minl1 && abs(Gy(edge-rows+l1)+1e-6)>=abs(Gy(edge-rows+l1-1)+1e-8)  
        l1 = l1 - 1;
    end
    while l2<maxl2 && abs(Gy(edge-rows+l2)+1e-6)>=abs(Gy(edge-rows+l2+1)+1e-8)  
        l2 = l2 + 1;
    end
    while m1>-4 && abs(Gy(edge+m1)+1e-6)>=abs(Gy(edge+m1-1)+1e-8)
        m1 = m1 - 1;
    end
    while m2<4 && abs(Gy(edge+m2)+1e-6)>=abs(Gy(edge+m2+1)+1e-8)
        m2 = m2 + 1;
    end
    while r1>minr1 && abs(Gy(edge+rows+r1)+1e-6)>=abs(Gy(edge+rows+r1-1)+1e-8)
        r1 = r1 - 1;
    end
    while r2<maxr2 && abs(Gy(edge+rows+r2)+1e-6)>=abs(Gy(edge+rows+r2+1)+1e-8)
        r2 = r2 + 1;
    end
    window = zeros(9,3);
    window(l1+5:l2+5,1) = 1;
    window(m1+5:m2+5,2) = 100;
    window(r1+5:r2+5,3) = 1;
    
    % compute intensities
    if m > 0
        AA = (G(edge+m2) + G(edge+rows+r2)) / 2;
        BB = (G(edge-rows+l1) + G(edge+m1)) / 2;
    else
        AA = (G(edge-rows+l2) + G(edge+m2)) / 2;
        BB = (G(edge+m1) + G(edge+rows+r1)) / 2;
    end

    % search for a second close edge
    uBorder = false;
    dBorder = false;
    if m1 > -4 
        partial = abs(G(edge+m1-2));
        if partial>abs(G(edge)/4) && partial>threshold/2
            uBorder = true;
        end
    end
    if m2 < 4 
        partial = abs(G(edge+m2+2));
        if partial>abs(G(edge)/4) && partial>threshold/2
            dBorder = true;
        end
    end
    SL=0; SM=0; SR=0;
    j = floor((edges(k)-1) / rows) + 1;
    i = edges(k) - rows*(j-1);

%     if i == 171 && j == 177
%         1 == 1;
%     end
% 
%     if k == 1575
%         1 == 1;
%     end

    if uBorder || dBorder
        rimvt = F(i-5:i+5, j-2:j+2);
        if uBorder
            if m > 0
                BB = (F(edge+m1) + F(edge-rows+l1)) / 2;
                p = 1;
            else
                BB = (F(edge+m1) + F(edge+rows+r1)) / 2;
                p = 0;
            end
            if Gy(edge-2*rows+l1+p) * Gy(edge) > 0
                ll = l1 + p - 1;
            else
                ll = l1 + p;
            end
            if Gy(edge+2*rows+r1+1-p) * Gy(edge) > 0
                rr = r1 - p;
            else
                rr = r1 + 1 - p;
            end
            rimvt(1:ll+6,1) = BB;
            rimvt(1:l1+6,2) = BB;
            rimvt(1:m1+6,3) = BB;
            rimvt(1:r1+6,4) = BB;
            rimvt(1:rr+6,5) = BB;
            l1=-3+m; m1=-3; r1=-3-m;
        end
        if dBorder
            if m > 0
                AA = (F(edge+m2) + F(edge+rows+r2)) / 2;
                p = 1;
            else
                AA = (F(edge+m2) + F(edge-rows+l2)) / 2;
                p = 0;
            end
            if Gy(edge-2*rows+l2+p-1) * Gy(edge) > 0
                ll = l2 + p;
            else
                ll = l2 + p - 1;
            end
            if Gy(edge+2*rows+r2-p) * Gy(edge) > 0
                rr = r2 + 1 - p;
            else
                rr = r2 - p;
            end
            rimvt(ll+6:11,1) = AA;
            rimvt(l2+6:11,2) = AA;
            rimvt(m2+6:11,3) = AA;
            rimvt(r2+6:11,4) = AA;
            rimvt(rr+6:11,5) = AA;
            l2=3+m; m2=3; r2=3-m;
        end
        rimv2 = (rimvt(1:9,1:3) + rimvt(1:9,2:4) + rimvt(1:9,3:5) + ...
            rimvt(2:10,1:3) + rimvt(2:10,2:4) + rimvt(2:10,3:5) + ...
            rimvt(3:11,1:3) + rimvt(3:11,2:4) + rimvt(3:11,3:5))/9;
        for n=l1+5:l2+5 SL = SL + rimv2(n,1); end;
        for n=m1+5:m2+5 SM = SM + rimv2(n,2); end;
        for n=r1+5:r2+5 SR = SR + rimv2(n,3); end;
    else  
        for n=l1:l2 SL = SL + G(edge-rows+n); end;
        for n=m1:m2 SM = SM + G(edge+n); end;
        for n=r1:r2 SR = SR + G(edge+rows+n); end;
    end
    
    % compute edge features
    den = 2 * (AA-BB);
    if order == 2 
        c(k) = (SL+SR-2*SM + AA*(2*m2-l2-r2) - BB*(2*m1-l1-r1)  + 1e-8) / den;
    else
        c(k) = 0; 
    end
    a(k) = (2*SM - AA*(1+2*m2) - BB*(1-2*m1)  + 1e-8) / den - w*c(k);
    if abs(a(k)) > 1
        valid(k) = false;
        continue;
    end
    valid(k) = true;
    b(k) = (SR-SL + AA*(l2-r2) - BB*(l1-r1)  + 1e-8) / den;
    A(k) = AA; 
    B(k) = BB;
    s = sign(AA-BB);
    nx(k) = s / sqrt(1+b(k)^2)*b(k);
    ny(k) = s / sqrt(1+b(k)^2);
    curv(k) = 2*c(k) / ((1+b(k)^2)^1.5);
    if Gy(edge)<0
        curv(k) = -curv(k);
    end
    
    % generate circle subimage
    R = abs(1/curv(k));
    if R > 1e4
        R = 1e4;
    end
    if R < 4.5
        R = 4.5;
    end
    if curv(k) > 0
        s = -1;
        innerIntensity = min(AA, BB);
        outerIntensity = max(AA, BB);
    else
        s = 1;
        innerIntensity = max(AA, BB);
        outerIntensity = min(AA, BB);
    end



    center = [x(edge)+s*R*nx(k), y(edge)-a(k)+s*R*ny(k)];
    subimage = circleVerticalWindow(j, i, center(1), center(2), ...
        R, innerIntensity, outerIntensity, pixelGridResol);
    
%     savename = "./subimage/subimage" + num2str(k) + ".txt";
%     fileID = fopen(savename,'w+');
%     [subimage_rows, subimage_cols] = size(subimage);
%     for i = 1:subimage_rows
%         for j = 1:subimage_cols
%             fprintf(fileID, '%.4f ', subimage(i,j));
%         end
%         fprintf(fileID, '\n');
%     end
%     fclose(fileID);

%     outwin = window.*subimage;
%     savename = "./outwin/outwin" + num2str(k) + ".txt";
%     fileID = fopen(savename,'w+');
%     [window_rows, window_cols] = size(outwin);
%     for iii = 1:window_rows
%         for jjj = 1:window_cols
%             fprintf(fileID, '%.4f ', outwin(iii,jjj));
%         end
%         fprintf(fileID, '\n');
%     end
%     fclose(fileID);

%     if i == 173 && j == 179
%         1 == 1;
%     end


    
    % update counter and intensity images
    I(i-4:i+4,j-1:j+1) = I(i-4:i+4,j-1:j+1)+window.*subimage;
    C(i-4:i+4,j-1:j+1) = C(i-4:i+4,j-1:j+1)+window;
end

% fileID = fopen('I.txt','w+');
% [window_rows, window_cols] = size(I);
% for iii = 1:window_rows
%     for jjj = 1:window_cols
%         fprintf(fileID, '%.4f ', I(iii,jjj));
%     end
%     fprintf(fileID, '\n');
% end
% fclose(fileID);


% fileID = fopen('valid.txt','w+');
% for iii = 1:size(valid,1)
%     fprintf(fileID, '%d\n', valid(iii));
% end
% fclose(fileID);

% remove invalid values
edges = edges(valid);
A = A(valid);
B = B(valid);

% save edges info
%n = ones(size(edges,1),1);
%n(Gy(edges)<0) = -1;
ep.position = [ep.position; edges];
ep.x = [ep.x; x(edges)];
ep.y = [ep.y; y(edges) - a(valid)];
ep.nx = [ep.nx; nx(valid)];
ep.ny = [ep.ny; ny(valid)];
ep.curv = [ep.curv; curv(valid)];
ep.i0 = [ep.i0; min(A,B)];
ep.i1 = [ep.i1; max(A,B)];


%% detecte edge pixels with maximum Gx (not including margins)
absGxInner = abs(Gx(3:rows-2,6:cols-5));
E = false(rows, cols);
E(3:rows-2,6:cols-5) = grad(3:rows-2,6:cols-5)>threshold & ...
    absGxInner>abs(Gy(3:rows-2,6:cols-5)+1e-8) & ...
    absGxInner>=abs(Gx(3:rows-2,5:cols-6)+1e-8) & ...
    absGxInner>abs(Gx(3:rows-2,7:cols-4)+1e-8);

% fileID = fopen('E.txt','w+');
% [window_rows, window_cols] = size(E);
% for iii = 1:window_rows
%     for jjj = 1:window_cols
%         fprintf(fileID, '%d ', E(iii,jjj));
%     end
%     fprintf(fileID, '\n');
% end
% fclose(fileID);

edges = (x(E)-1)*rows+y(E);
A = zeros(size(edges,1),1);
B = zeros(size(edges,1),1);
a = zeros(size(edges,1),1);
b = zeros(size(edges,1),1);
c = zeros(size(edges,1),1);
nx = zeros(size(edges,1),1);
ny = zeros(size(edges,1),1);
curv = zeros(size(edges,1),1);
valid = logical(zeros(size(edges,1),1));

% compute all horizontal edges
for k=1:size(edges,1)
    edge = edges(k);
 
    % compute window floating limits
    m1=-1;  m2=1;
    if Gx(edge)*Gy(edge)>=0
        m = 1;
        l1=-1;      r2=1;
        minl1=-3;   maxr2=3;
        l2=1;       r1=-1;
        maxl2=4;    minr1=-4;
    else
        m = -1;
        l1=-1;      r2=1;
        minl1=-4;   maxr2=4;
        l2=1;       r1=-1;
        maxl2=3;    minr1=-3;
    end        
    while l1>minl1 && abs(Gx(edge-1+l1*rows)+1e-6)>=abs(Gx(edge-1+(l1-1)*rows)+1e-8)  
        l1 = l1 - 1;
    end
    while l2<maxl2 && abs(Gx(edge-1+l2*rows)+1e-6)>=abs(Gx(edge-1+(l2+1)*rows)+1e-8)  
        l2 = l2 + 1;
    end
    while m1>-4 && abs(Gx(edge+m1*rows)+1e-6)>=abs(Gx(edge+(m1-1)*rows)+1e-8)
        m1 = m1 - 1;
    end
    while m2<4 && abs(Gx(edge+m2*rows)+1e-6)>=abs(Gx(edge+(m2+1)*rows)+1e-8)
        m2 = m2 + 1;
    end
    while r1>minr1 && abs(Gx(edge+1+r1*rows)+1e-6)>=abs(Gx(edge+1+(r1-1)*rows)+1e-8)
        r1 = r1 - 1;
    end
    while r2<maxr2 && abs(Gx(edge+1+r2*rows)+1e-6)>=abs(Gx(edge+1+(r2+1)*rows)+1e-8)
        r2 = r2 + 1;
    end
    window = zeros(3,9);
    window(1,l1+5:l2+5) = 1;
    window(2,m1+5:m2+5) = 100;
    window(3,r1+5:r2+5) = 1;

    % compute intensities
    if m > 0
        AA = (G(edge+m2*rows) + G(edge+1+r2*rows)) / 2;
        BB = (G(edge-1+l1*rows) + G(edge+m1*rows)) / 2;
    else
        AA = (G(edge-1+l2*rows) + G(edge+m2*rows)) / 2;
        BB = (G(edge+m1*rows) + G(edge+1+r1*rows)) / 2;
    end

    % search for a second close edge
    uBorder = false;
    dBorder = false;
    if m1 > -4 
        partial = abs(G(edge+(m1-2)*rows));
        if partial>abs(G(edge)/4) && partial>threshold/2
            uBorder = true;
        end
    end
    if m2 < 4 
        partial = abs(G(edge+(m2+2)*rows));
        if partial>abs(G(edge)/4) && partial>threshold/2
            dBorder = true;
        end
    end
    SL=0; SM=0; SR=0;
    j = floor((edges(k)-1) / rows) + 1;
    i = edges(k) - rows*(j-1);

    if (k == 8)
        1 == 1;
    end

    if uBorder || dBorder
        rimvt = F(i-2:i+2, j-5:j+5);
        if uBorder
            if m > 0
                BB = (F(edge+m1*rows) + F(edge-1+l1*rows)) / 2;
                p = 1;
            else
                BB = (F(edge+m1*rows) + F(edge+1+r1*rows)) / 2;
                p = 0;
            end
            if Gx(edge-2+(l1+p)*rows) * Gx(edge) > 0
                ll = l1 + p - 1;
            else
                ll = l1 + p;
            end
            if Gx(edge+2+(r1+1-p)*rows) * Gx(edge) > 0
                rr = r1 - p;
            else
                rr = r1 + 1 - p;
            end
            rimvt(1,1:ll+6) = BB;
            rimvt(2,1:l1+6) = BB;
            rimvt(3,1:m1+6) = BB;
            rimvt(4,1:r1+6) = BB;
            rimvt(5,1:rr+6) = BB;
            l1=-3+m; m1=-3; r1=-3-m;
        end
        if dBorder
            if m > 0
                AA = (F(edge+m2*rows) + F(edge+1+r2*rows)) / 2;
                p = 1;
            else
                AA = (F(edge+m2*rows) + F(edge-1+l2*rows)) / 2;
                p = 0;
            end
            if Gx(edge-2+(l2+p-1)*rows) * Gx(edge) > 0
                ll = l2 + p;
            else
                ll = l2 + p - 1;
            end
            if Gx(edge+2+(r2-p)*rows) * Gx(edge) > 0
                rr = r2 + 1 - p;
            else
                rr = r2 - p;
            end
            rimvt(1,ll+6:11) = AA;
            rimvt(2,l2+6:11) = AA;
            rimvt(3,m2+6:11) = AA;
            rimvt(4,r2+6:11) = AA;
            rimvt(5,rr+6:11) = AA;
            l2=3+m; m2=3; r2=3-m;
        end
        rimv2 = (rimvt(1:3,1:9) + rimvt(2:4,1:9) + rimvt(3:5,1:9) + ...
            rimvt(1:3,2:10) + rimvt(2:4,2:10) + rimvt(3:5,2:10) + ...
            rimvt(1:3,3:11) + rimvt(2:4,3:11) + rimvt(3:5,3:11))/9;
        for n=l1+5:l2+5 SL = SL + rimv2(1,n); end;
        for n=m1+5:m2+5 SM = SM + rimv2(2,n); end;
        for n=r1+5:r2+5 SR = SR + rimv2(3,n); end;
    else  
        for n=l1:l2 SL = SL + G(edge-1+n*rows); end;
        for n=m1:m2 SM = SM + G(edge+n*rows); end;
        for n=r1:r2 SR = SR + G(edge+1+n*rows); end;
    end
  
    % compute edge features
    den = 2 * (AA-BB);
    if order == 2
        c(k) = (SL+SR-2*SM + AA*(2*m2-l2-r2) - BB*(2*m1-l1-r1)+ 1e-8) / den;
    else
        c(k) = 0; 
    end
    a(k) = (2*SM - AA*(1+2*m2) - BB*(1-2*m1)+ 1e-8) / den - w*c(k);
    if abs(a(k)) > 1
        valid(k) = false;
        continue;
    end
    valid(k) = true;
    b(k) = (SR-SL + AA*(l2-r2) - BB*(l1-r1)+ 1e-8) / den;
    A(k) = AA; 
    B(k) = BB;
    s = sign(AA-BB);
    nx(k) = s / sqrt(1+b(k)^2);
    ny(k) = s / sqrt(1+b(k)^2)*b(k);
    curv(k) = 2*c(k) / ((1+b(k)^2)^1.5);
    if Gx(edge)<0
        curv(k) = -curv(k);
    end
    
    % generate circle subimage
    R = abs(1/curv(k));
    if R > 1e4
        R = 1e4;
    end
    if R < 4.5
        R = 4.5;
    end
    if curv(k) > 0
        s = -1;
        innerIntensity = min(AA, BB);
        outerIntensity = max(AA, BB);
    else
        s = 1;
        innerIntensity = max(AA, BB);
        outerIntensity = min(AA, BB);
    end
    center = [x(edge)-a(k)+s*R*nx(k), y(edge)+s*R*ny(k)];
    subimage = circleHorizontalWindow(j, i, center(1), center(2), ...
        R, innerIntensity, outerIntensity, pixelGridResol);
    
    % update counter and intensity images
    I(i-1:i+1,j-4:j+4) = I(i-1:i+1,j-4:j+4)+window.*subimage;
    C(i-1:i+1,j-4:j+4) = C(i-1:i+1,j-4:j+4)+window;
end




% remove invalid values
edges = edges(valid);
A = A(valid);
B = B(valid);

% save edges info
ep.position = [ep.position; edges];
ep.x = [ep.x; x(edges) - a(valid)];
ep.y = [ep.y; y(edges)];
ep.nx = [ep.nx; nx(valid)];
ep.ny = [ep.ny; ny(valid)];
ep.curv = [ep.curv; curv(valid)];
ep.i0 = [ep.i0; min(A,B)];
ep.i1 = [ep.i1; max(A,B)];

mkdir('./results/');
fileID = fopen('./results/position.txt','w+');
for iii = 1:size(ep.position,1)
    fprintf(fileID, '%d\n', ep.position(iii));
end
fclose(fileID);

fileID = fopen('./results/x.txt','w+');
for iii = 1:size(ep.x,1)
    fprintf(fileID, '%.4f\n', ep.x(iii));
end
fclose(fileID);

fileID = fopen('./results/y.txt','w+');
for iii = 1:size(ep.y,1)
    fprintf(fileID, '%.4f\n', ep.y(iii));
end
fclose(fileID);

fileID = fopen('./results/nx.txt','w+');
for iii = 1:size(ep.nx,1)
    fprintf(fileID, '%.4f\n', ep.nx(iii));
end
fclose(fileID);

fileID = fopen('./results/ny.txt','w+');
for iii = 1:size(ep.ny,1)
    fprintf(fileID, '%.4f\n', ep.ny(iii));
end
fclose(fileID);

fileID = fopen('./results/curv.txt','w+');
for iii = 1:size(ep.curv,1)
    fprintf(fileID, '%.4f\n', ep.curv(iii));
end
fclose(fileID);

fileID = fopen('./results/i0.txt','w+');
for iii = 1:size(ep.i0,1)
    fprintf(fileID, '%.4f\n', ep.i0(iii));
end
fclose(fileID);

fileID = fopen('./results/i1.txt','w+');
for iii = 1:size(ep.i1,1)
    fprintf(fileID, '%.4f\n', ep.i1(iii));
end
fclose(fileID);


%% compute final subimage
I1 = I(C>0);
C1 = C(C>0);
I(C>0) = I(C>0) ./ C(C>0);
I(C==0) = G(C==0);

fileID = fopen('./results/I.txt','w+');
[window_rows, window_cols] = size(I);
for iii = 1:window_rows
    for jjj = 1:window_cols
        fprintf(fileID, '%.4f ', I(iii,jjj));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

end

