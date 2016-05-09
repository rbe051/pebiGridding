function [p] = clippedPebi2D(p, bnd)

dt = delaunayTriangulation(p);
m = size(bnd,1);
dtB.Points = bnd;
dtB.ConnectivityList = reshape([1:m-1;2:m],[],1);

[V, C, symV] = clipGrid(dt, dtB);


end