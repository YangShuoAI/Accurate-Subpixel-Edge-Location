function G = smooth(F)

% smooth image
G = double(F);
[rows, cols] = size(F);
G(2:rows-1,2:cols-1) = ...
    (F(1:rows-2,1:cols-2) + F(1:rows-2,2:cols-1) + F(1:rows-2,3:cols) + ...
    F(2:rows-1,1:cols-2) + F(2:rows-1,2:cols-1) + F(2:rows-1,3:cols) + ...
    F(3:rows,1:cols-2) + F(3:rows,2:cols-1) + F(3:rows,3:cols))/9;

% show image
imshow(uint8(G));
axis('image');
end