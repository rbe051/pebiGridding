function [newPoints, dt] = eqInterpret(path, dt)
    % Interpolate a path with equiv distant points
    % Arguments:
    %   path       n*2 array of coordinates of points which are to be
    %              interpolated
    %   dt         distance between interpolation points
    %
    % Return:
    %   newPoints  The interpolated points
    %   dt         distance between the interpolated points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Runar Lie Berge (runarlb@stud.ntnu.no)
    % January 2016
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    linesDist = sqrt(sum(diff(path,[],1).^2,2));
    linesDist = [0; linesDist]; % add the starting point
    cumDist = cumsum(linesDist);
    dt = cumDist(end)/ceil(cumDist(end)/dt);
    newPointsLoc = 0:dt:cumDist(end);
        
    newPoints = interp1(cumDist, path, newPointsLoc);    
end
