function [res] = mlqt(pts, bndr, gridSize, level, levelTol, distTol)
    n = size(bndr,1);
    repPnt = repmat(pts,n,1);
    if any(sum((repPnt-bndr).^2,2) < distTol^2) && level<= levelTol
        shift = gridSize/2;
        res = [mlqt(pts + [shift,shift], bndr, gridSize/2, level+1, levelTol, distTol/2);...
                  mlqt(pts + [shift,-shift], bndr, gridSize/2, level+1, levelTol, distTol/2);...
                  mlqt(pts + [-shift,-shift], bndr, gridSize/2, level+1, levelTol, distTol/2);...
                  mlqt(pts + [-shift,shift], bndr, gridSize/2, level+1, levelTol, distTol/2)];
            
    else
        res = {pts, gridSize*(1-10^-6)};
    end
end