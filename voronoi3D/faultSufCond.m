function [pts,removed] = faultSufCond(pts, CC, CR)
    TOL = 10*eps;

    nc = size(CC,1);
    np = size(pts,1);
    
    CRSqr = CR.^2;
    removed = zeros(np,1);
    for i = 1:nc
      distSqr = sum(bsxfun(@minus, CC(i,:), pts).^2,2);
      removed = removed + (distSqr<CRSqr(i)-TOL);
    end
    removed = logical(removed);
    pts = pts(~removed,:);
end
