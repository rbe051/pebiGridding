function [] = plotWells(G)

if length(strmatch('computeGeometry',G.type,'exact'))<1
    G = computeGeometry(G);
end

centroids = G.cells.centroids(G.cells.isWell,:);
plot(centroids(:,1),centroids(:,2), '.');
end
