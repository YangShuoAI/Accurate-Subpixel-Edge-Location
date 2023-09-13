function visEdges(edges, varargin)
%VISEDGES displays a list of edges
%
%   VISEDGES(EDGES) displays a list of EDGES in the same format as the 
%   result of subpixelEdges method.
%
%   VISEDGES(EDGES,PARAM1,VAL1,PARAM2,VAL2,...) displays the list of EDGES
%   using name-value pairs to control what must be visualized. These 
%   parameters include:
%
%   'showEdges' - Specifies if edges are displayed:
%               true:  edges are displayed (default)
%               false:  edges are not displayed 
% 
%   'showNormals' - Specifies if normal vectors ared displayed:
%               true:  normal vectors are displayed (default)
%               false:  normal vectors are not displayed 
%
%   'condition' - Specifies what edges to draw. Must be a logical vector of
%   the same size of edges

hold on;
showEdges = true;
showNormals = true;
condition = [];

%% parse optional input parameters
v = 1;
while v < numel(varargin)
    switch varargin{v}
        case 'showEdges'
            assert(v+1<=numel(varargin));
            showEdges = varargin{v+1};
        case 'showNormals'
            assert(v+1<=numel(varargin));
            showNormals = varargin{v+1};
        case 'condition'
            assert(v+1<=numel(varargin));
            condition = varargin{v+1};
            if length(condition) ~= length(edges.position)
                error('condition must be the same size than edges');
            end
        otherwise
            error('Unsupported parameter: %s',varargin{v});
    end
    v = v+2;
end

%% apply condition if needed
if ~isempty(condition)
    edges.position = edges.position(condition);
    edges.x = edges.x(condition);
    edges.y = edges.y(condition);
    edges.nx = edges.nx(condition);
    edges.ny = edges.ny(condition);
    edges.curv = edges.curv(condition);
    edges.i0 = edges.i0(condition);
    edges.i1 = edges.i1(condition);
end

%% display edges
if showEdges
    seg = 0.6;
    quiver(edges.x-seg/2*edges.ny, edges.y+seg/2*edges.nx, ...
        seg*edges.ny, -seg*edges.nx, 0, 'r.');
end

%% display normal vectors
if showNormals
    quiver(edges.x, edges.y, edges.nx, edges.ny, 0, 'b');
end

hold off
