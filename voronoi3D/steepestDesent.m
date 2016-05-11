function [x, f, gNorm] = steepestDesent(x0, F, dt, varargin)
    % limitet-memory bfgs optimization function
    %
    % Arguments:
    %   x0          initial guess
    %   F           Objective function
    %   dt          delaunay triangulation class of the boundary
    %
    % varargin:
    %   storedVec   Number of vectors used to approximate hessian
    %   maxIt       Maximum number of iterations
    %   tol         Convergence tolerance. Convergence test is
    %               gradF(x) <= tol*gradF(x0)
    %   minStep     if steplength is less than minStep, the function
    %               returns
    % Returns:
    %   x           optimal point
    %   f           functioin value at each step
    %   gNorm       norm of the gradient at each step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Runar Lie Berge (runarlb@stud.ntnu.no)                           2016
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


opt = struct('maxIt',     1000, ...
             'tol',       1e-6, ...
             'minStep',   10*eps,...
             'betaM',      0.5, ...
             'betaP',      1.1,...
             's0',           0.1);

opt    = merge_options(opt,varargin{:});

maxIt  = opt.maxIt;
tol    = opt.tol;
betaM  = opt.betaM;
betaP  = opt.betaP;

x      = x0;
s      = opt.s0;
[f,g]  = F(x0);
gNorm  = norm(g,2);

for k = 1:maxIt
  % Find search direction
  [f(k+1), gNew] = F(x);
  gNorm(k+1) = norm(g,2);
  if gNorm(k+1) <= tol*gNorm(1)
    return
  end
  p = -gNew/norm(gNorm(k+1));
  if dot(g,gNew)<0
    s = s*betaM;
  else
    s = s*betaP;
  end

  x0 = x;
  x = x + s*p;
  g = gNew;
  fprintf('%3d: f=%10.3e, |df|=%10.3e, |xk+1-xk| = %10.3e\n', ...
           k, f(k+1), gNorm(k+1), norm(x-x0));
end