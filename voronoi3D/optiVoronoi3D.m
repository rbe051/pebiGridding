function [G,optPts,f,g] = optiVoronoi3D(pts, bndr, varargin)
    % Creates an optimal voronoi voronoi diagram. By optimal means that all
    % voronoi seeds coincide with the cell mass centers.
    %
    % Arguments:
    %   pts         n x 3 array of initial of seeds points
    %   bndr        nb x 3 array of points defining the convext boundary
    % Varagin:
    %   density     density function in the domain
    %   storedVec   number of stored vectors used to approximate the
    %               hessian in the lbfgs algorithm
    %   tol         convergence tolerance
    %   maxIt       maximum number of iteration in the lbfgs algorithm
    %   minStep     if step length is smaler than minStep, the lbfgs
    %               algorithm stops
    % Returns:
    %   G           mrst grid structure
    %   optPts      optimal grid seeds
    %   f           function value at each step
    %   g           gradient value at each step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Runar Lie Berge                                                  2016
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    assert (size(pts, 2) == 3, ...
          ['Function ''%s'' is only supported in three ', ...
           'space dimensions.'], mfilename);
    assert (size(bndr, 2) == 3, ...
          ['Function ''%s'' is only supported in three ', ...
           'space dimensions.'], mfilename);

    opt = struct('density',   @(x) ones(size(x,1),1),...
                 'storedVec', 5,                     ...
                 'tol',       1e-6,                 ...         
                 'maxIt',     1000,                  ...
                 'minStep',   10*eps);
    opt = merge_options(opt, varargin{:});
    
    F = @(pts) objectiveFunc(pts, bndr, opt.density);

    pts = reshape(pts',[],1);
    dt = delaunayTriangulation(bndr);
    [optPts, f, g] = lbfgs(pts, F, dt, 'storedVec',opt.storedVec, ...
                                       'maxIt',    opt.maxIt,...
                                       'tol',      opt.tol);

    optPts = reshape(optPts,3,[])';
    G = voronoi3D(optPts, bndr);
end

function [f, g] = objectiveFunc(pts, bndr, rho)
    pts = reshape(pts,3,[])';

    G = voronoi3D(pts, bndr);
    G = computeGeometry(G);

    intFun = @(x,i) sum(repmat(rho(x),1,3).*(x-repmat(pts(i,:),size(x,1),1)).^2,2);
    f = sum(polyhedronInt(G,1:G.cells.num, intFun));

    massFun = @(x,i) rho(x);
    masses = polyhedronInt(G,1:G.cells.num,massFun);
    g = reshape((2*repmat(masses,1,3).*(pts - G.cells.centroids))',[],1);

end
