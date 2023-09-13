classdef EdgePixel    
    properties
        position;       % 1D index inside image 
        x, y;           % subpixel position
        nx, ny;         % normal vector (normalized)
        curv;           % curvature
        i0, i1;         % intensities at both sides
    end
    
    methods
    end 
end

