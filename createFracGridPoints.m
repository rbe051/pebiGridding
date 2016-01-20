function [Pts, gridSpacing, circCenter, circRadius, CCid] = ...
    createFracGridPoints(faultLine, fracDs, circleFactor, varargin) 
    
    opt = struct('distFunc', @huniform);

    opt = merge_options(opt,varargin{:});
    fh = opt.distFunc;
    assert(0.5<circleFactor && circleFactor<1)
    assert(size(faultLine,1)>1 && size(faultLine,2)==2);

    %[circCenter, ~] = eqInterpret(faultLine, fracDs, h);
    circCenter = interFaultLine(faultLine, fh, fracDs);
    %This is equidistante if you follow the line described by fracLine,
    %but the new points may be a bit to close if the fracture has
    %sharp corners and/or you upsample the line.

    numOfFracPts = size(circCenter,1)-1;
    if numOfFracPts <= 0
        Pts = [];
        gridSpacing = [];
        return
    end

    % Create left vector and right vector(to contain fracture points).
    left = zeros(numOfFracPts, 2);
    right = zeros(numOfFracPts, 2);
    ptGrSp = zeros(numOfFracPts, 1);
    lineLength = sqrt(sum((circCenter(2:end,:)-circCenter(1:end-1,:)).^2, 2));
    
    for j=1:numOfFracPts
        %lineLength = norm(fracLine(j+1,:) - fracLine(j,:));
        %Because the line lengths are not completely uniform we have to
        %calculate this for each segment.
        fractureRadius = lineLength(j)/2*sqrt(4*circleFactor^2 -1);
        n1 = (circCenter(j+1,:)-circCenter(j,:))/lineLength(j);   %Unit vector
        n2 = [-n1(2), n1(1)];                              %Unit normal
        left(j,:) = circCenter(j,:) + lineLength(j)/2*n1 + fractureRadius*n2;
        right(j,:) = circCenter(j,:) + lineLength(j)/2*n1 - fractureRadius*n2;
        ptGrSp(j) = 2*fractureRadius;
    end
    circRadius = circleFactor*[lineLength;lineLength(end)];
    Pts = [right;left];
    CCid = [1:size(left,1),1:size(right,1)]';
    gridSpacing = [ptGrSp; ptGrSp];
end


function [newPoints, dt] = eqInterpret(path, dt)
    linesDist = sqrt(sum(diff(path,[],1).^2,2));
    linesDist = [0; linesDist]; % add the starting point
    cumDist = cumsum(linesDist);
    dt = cumDist(end)/ceil(cumDist(end)/dt);
    newPointsLoc = 0:dt:cumDist(end);
    if length(newPointsLoc)<2 % Fracture to short
        error('Fracture grid size larger than fracture')
    end
        
    newPoints = interp1(cumDist, path, newPointsLoc);
      
end

