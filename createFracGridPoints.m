function [Pts, gridSpacing] = createFracGridPoints(fracLine, fracDs, circleFactor) 
    assert(0.5<circleFactor && circleFactor < 1)
    assert(0<fracDs)
    assert(size(fracLine,1)>1, size(fracLine,2)==2);

    [fracLine, fracDs] = eqInterpret(fracLine, fracDs);
    %This is equidistante if you follow the line described by fracLine,
    %but the new points may be a bit to close if the fracture has
    %sharp corners and/or you upsample the line.

    numOfFracPts = size(fracLine,1)-1;

    % Create left vector and right vector(to contain fracture points).
    left = zeros(numOfFracPts, 2);
    right = zeros(numOfFracPts, 2);
    ptGrSp = zeros(numOfFracPts, 1);
    
    for j=1:numOfFracPts
        lineLength = norm(fracLine(j+1,:) - fracLine(j,:));
        %Because the line lengths are not completely uniform we have to
        %calculate this for each segment.
        fractureRadius = lineLength/2*sqrt(4*circleFactor^2 -1);
        n1 = (fracLine(j+1,:)-fracLine(j,:))/lineLength;   %Unit vector
        n2 = [-n1(2), n1(1)];                              %Unit normal
        left(j,:) = fracLine(j,:) + lineLength/2*n1 + fractureRadius*n2;
        right(j,:) = fracLine(j,:) + lineLength/2*n1 - fractureRadius*n2;
        ptGrSp(j) = (2-10^-6*fracDs)*fractureRadius;
    end
    
    Pts = [right;left];
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




