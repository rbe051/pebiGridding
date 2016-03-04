function [Pts, gridSpacing] = createWellGridPoints(wellLine, wellDs) 
    % Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
    assert(0<wellDs)
    assert(size(wellLine,1)>=1, size(wellLine,2)==2);
    if (size(wellLine,1) == 1)
        Pts = wellLine;
        gridSpacing = wellDs*(1-10^-6);
    else
        [Pts, ~] = eqInterpret(wellLine, wellDs);
        gridSpacing = sqrt(sum(diff(Pts,1,1).^2,2));
        gridSpacing = [gridSpacing(1);gridSpacing;gridSpacing(end)];
        gridSpacing = min([gridSpacing(1:end-1), gridSpacing(2:end)],[],2);
    end
end
