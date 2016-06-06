function [G,optPts,f,g] = createCVD(pts,bnd,varargin)

opt = struct('fixedPts', [],...
             'rho',     @(pts) ones(size(pts,1),1));
opt = merge_options(opt,varargin{:});
fixedPts = opt.fixedPts;
rho = opt.rho;

nf = size(fixedPts,1);

dt = delaunayTriangulation(bnd);
p = [fixedPts;pts];
F = @(pts) objFunc(pts, bnd,nf,rho);

p = reshape(p',[],1);
[optPts,f,g] = lbfgs(p, F, dt,'tol',1e-5,'storedVec',7);
optPts = reshape(optPts,2,[])';

G = clippedPebi2D(optPts, bnd);

end


function [f, g] = objFunc(p, bnd,nf,rho)
    pts = reshape(p,2,[])';

    G = clippedPebi2D(pts, bnd);
    %G  = mirrorPebi2D(pts, bnd);
    G = computeGeometry(G);
    G = mrstGridWithFullMappings(G);


    intFun = @(x,i) sum(repmat(rho(x),1,2).*(x-repmat(pts(i,:),size(x,1),1)).^2,2);
    f = sum(polygonInt_v2(G,1:G.cells.num, intFun,7));

    massFun = @(x,i) rho(x);
    masses = polygonInt_v2(G,1:G.cells.num,massFun,7);
    g = reshape((2*repmat(masses,1,2).*(pts - G.cells.centroids))',[],1);
    g(1:2*nf) = 0;
end