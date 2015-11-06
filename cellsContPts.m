function [c] = cellsContPts(G, pts)

if length(strmatch('computeGeometry',G.type,'exact'))<1
    G = computeGeometry(G);
end

n = G.cells.num;
c = zeros(n,1);
for i = 1:n
    for j = 1:size(pts,1)
        faces = G.cells.faces(G.cells.facePos(i):(G.cells.facePos(i+1)-1));
        normals = G.faces.normals(faces,:);
        % Normals should point out of cell:
        normals = normals.*repmat((-1+2*(G.faces.neighbors(faces',1)==i)),1,2);
        centroids = G.faces.centroids(faces,:);

        c(i)=c(i) + all(diag((repmat(pts(j,:),size(centroids,1),1)-centroids)*normals'<=0));
    end
end
c = c>0;
end