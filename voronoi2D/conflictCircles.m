function [I] = conflictCircles(Pts, CC, CR)
    TOL = 10*eps;

    nc = size(CC,1);
    np = size(Pts,1);
    removed = zeros(np,1);
    CRSqr = CR.^2;
    I = cell(numel(CR),1);
    for i = 1:nc
        distSqr = sum((repmat(CC(i,:),np,1)-Pts).^2,2);
        I{i} = find(distSqr<CRSqr(i)-10*eps);
    end                
end
