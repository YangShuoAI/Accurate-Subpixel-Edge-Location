function [edges, RI] = subpixelEdges(image, threshold, varargin)
%SUBPIXELEDGES computes a subpixel edge detection using method published
%   in the following paper: 
%   "Accurate Subpixel Edge Location Based on Partial Area Effect"
%   http://www.sciencedirect.com/science/article/pii/S0262885612001850
%
%   EDGES = SUBPIXELEDGES(IMAGE,THRESHOLD) detects subpixel features for
%   each pixel belonging to an edge in IMAGE. IMAGE must be a grayscale
%   image. THRESHOLD specifies the minimum difference of intensity at both
%   sides of a pixel to be considered as an edge.
%
%   EDGES = SUBPIXELEDGES(IMAGE,THRESHOLD,PARAM1,VAL1,PARAM2,VAL2,...)
%   applies the edge detection using name-value pairs to control aspects of
%   the method. These parameters include:
%
%   'Order' - Specifies the order of the edge to find:
%               1:  firt order edges (straigh lines) 
%               2:  second order edges (default)
% 
%   'SmoothingIter' - Specifies how many smoothing iterations are needed to
%                       find final edges:
%               0:  Oriented to noise free images. No previous smoothing on
%                   the image. The detection is applied on the original 
%                   image values (section 3 of the paper).
%               1:  Oriented to low-noise images. The detection is applied 
%                   on the image previously smoothed by a 3x3 mask 
%                   (default) (sections 4 and 5 of the paper)
%               >1: Oriented to high-noise images. Several stages of 
%                   smoothing + detection + sinthetic image creation are 
%                   applied (section 6 of the paper). A few iterations are
%                   normally enough.
%
%   [EDGES, RI] = SUBPIXELEDGES(...) returns also the restored image RI
%   obtained in the last iteration (when iterations > 0) (see example5.m)
%
%   To see the intermediate results of every iteration, include every
%   single iteration inside a loop in the following way:
%       RI = image;
%       for n=1:iter
%           [edges, RI] = subpixelEdges(RI, threshold); 
%       end
%
%   EDGES is an array of structures with length equal to the number of
%   detected edge pixels in IMAGE. The fields of the structure denote
%   different features for each pixel:
%       position1D:     % index inside image 
%       x, y:           % subpixel position
%       nx, ny:         % normal vector (normalized)
%       curv:           % curvature
%       i0, i1:         % intensities at both sides
%

%% defaults for optional parameters
order = 2;
smoothingIter = 3;
if nargin<2
    error('Minimum two parameters needed (image and threshold)');
end

%% parse optional input parameters
v = 1;
while v < numel(varargin)
    switch varargin{v}
        case 'Order'
            assert(v+1<=numel(varargin));
            order = varargin{v+1};
            if (order<1 || order>2)
                error('Order must be 1 or 2');
            end
        case 'SmoothingIter'
            assert(v+1<=numel(varargin));
            smoothingIter = varargin{v+1};
        otherwise
            error('Unsupported parameter: %s',varargin{v});
    end
    v = v+2;
end

%% smooth is required when restored image is requested
if nargout==2 && smoothingIter<1
    error('smoothing is required if restored image is requested');
end

%% convert image to double
if strcmp(class(image), 'double') == false
    image = 255 * im2double(image);
end

%% call appropiated method
%addpath('EdgeDetector');
if smoothingIter==0
    edges = finalDetectorIter0(image, threshold, order);
elseif smoothingIter==1 && nargout==1
    edges = finalDetectorIter1(image, threshold, order);
else
    addpath('Synthetic');
    RI = image;
    for n=1:smoothingIter
        fprintf ('Iteration %d / %d...\n', n, smoothingIter);
        [edges, RI] = finalDetectorIterN(RI, threshold, order);
    end
end






