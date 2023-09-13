function edgesC = subsetEdges(edges, condition)
%SUBSETEDGES extracts a subset of edges that fulfill condition
%
%   Example:
%   edgesC = subsetEdges(edges, edges.x<10 & edges.y>30);

edgesC = EdgePixel; 
edgesC.position = edges.position(condition);
edgesC.x = edges.x(condition); 
edgesC.y = edges.y(condition); 
edgesC.nx = edges.nx(condition); 
edgesC.ny = edges.ny(condition); 
edgesC.curv = edges.curv(condition); 
edgesC.i0 = edges.i0(condition); 
edgesC.i1 = edges.i1(condition); 
