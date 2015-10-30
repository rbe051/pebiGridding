function [newPoints, dt] = eqInterpret(path, dt)
    linesDist = sqrt(sum(diff(path,[],1).^2,2));
    linesDist = [0; linesDist]; % add the starting point
    cumDist = cumsum(linesDist);
    dt = cumDist(end)/ceil(cumDist(end)/dt);
    newPointsLoc = 0:dt:cumDist(end);
        
    newPoints = interp1(cumDist, path, newPointsLoc);    
end
