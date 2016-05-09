function [G,optPts,f,g] = createCVD(k, n, x0)

 dt = delaunayTriangulation(x0);
% row1 = [linspace(0.165,0.835,3)',repmat(0.165,3,1)];
% row2 = [linspace(.12,.88,4)',repmat(0.5,4,1)];
% row3 = [linspace(0.165,0.835,3)',repmat(0.835,3,1)];
% pts = [row1;row2;row3];
pts = rand(k,2);

F = @(pts) objFunc(pts, n, x0);

pts = reshape(pts',[],1);
[optPts,f,g] = lbfgs(pts, F, dt,'tol',1e-7,'storedVec',10);
optPts = reshape(optPts,2,[])';


G = clippedPebi2D(optPts, n,x0);

end


function [f, g] = objFunc(p, n, x0)
    pts = reshape(p,2,[])';

    G = clippedPebi2D(pts, n, x0);
    G = computeGeometry(G);
    G = mrstGridWithFullMappings(G);

    rho = @(x) ones(size(x,1),1);
    intFun = @(x,i) sum(repmat(rho(x),1,2).*(x-repmat(pts(i,:),size(x,1),1)).^2,2);
    f = sum(polygonInt_v2(G,1:G.cells.num, intFun,7));

    massFun = @(x,i) ones(size(x,1),1);
    masses = polygonInt_v2(G,1:G.cells.num,massFun,7);
    g = reshape((2*repmat(masses,1,2).*(pts - G.cells.centroids))',[],1);

end