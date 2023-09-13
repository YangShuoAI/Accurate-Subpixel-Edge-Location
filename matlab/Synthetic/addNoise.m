function N = addNoise (I, m)
%NOISE add white gaussian noise of magnitude m to image I
% example m=0.20 --> noise of 20%
% inspired by
% http://stackoverflow.com/questions/16008228/
% using-imnoise-to-add-gaussian-noise-to-an-image

i = I / 255;
% n = imnoise(i, 'gaussian', 0, m*var(i(:)));
n = imnoise(i, 'gaussian', 0, 0.001*m);
N = n * 255;

% show image
%imshow(uint8(N));
%axis('image');
end