function [G,f,g] = cvdResPtsObjectiveFunc(p,bnd)


F = @(p) objFunc(p,bnd);

p = reshape(p',[],1);

[optP,f,g] = lbfgs(p, F, dt,'tol',1e-2,'storedVec',10);
optP = reshape(optP,2,[])';

G = clippedPebi2D(optP,bnd);
end


