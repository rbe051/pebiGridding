function [res] = mlqt(cellCenter, bndr, cellSize, varargin)
    
    opt = struct('level', 1, ...
                 'maxLev', 0, ...
                 'distTol', -1);
             
    opt = merge_options(opt, varargin{:});
    level = opt.level;
    maxLev = opt.maxLev;
    
    if level> maxLev
        res = {cellCenter, cellSize*(1-10^-6)};
        return
    end
    
    assert(cellSize>0);
    assert(size(opt.distTol,2) ==1 && size(opt.distTol,1) >0)
    
    if size(opt.distTol,1)==1
        if opt.distTol <= 0
            distTol = 1.5*cellSize/2;
        else
            distTol = opt.distTol;
        end
        distNext = distTol/2;
    else
        distTol = opt.distTol(level);
        distNext = opt.distTol;
    end
    
    n = size(bndr,1);
    repPnt = repmat(cellCenter,n,1);
    if any(max(abs(repPnt-bndr),[],2) < distTol)
        shift = cellSize/4;
        varArg = {'level', level+1, 'maxLev', maxLev, 'distTol', distNext};
        res = [mlqt(cellCenter + [shift,shift], bndr, cellSize/2, varArg{:});...
               mlqt(cellCenter + [shift,-shift], bndr, cellSize/2, varArg{:});...
               mlqt(cellCenter + [-shift,-shift], bndr, cellSize/2, varArg{:});...
               mlqt(cellCenter + [-shift,shift], bndr, cellSize/2, varArg{:})];            
    else
        res = {cellCenter, cellSize*(1-10^-6)};
    end
end