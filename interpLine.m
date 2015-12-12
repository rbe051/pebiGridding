function [newPoints] = interpLine(path, dt)
    distS = sqrt(sum(diff(path).^2,2));
    t = [0; cumsum(distS)];
    
    newPtsEval = [0; cumsum(dt)];
    onPath = newPtsEval < t(end); % Find points on path
    newPtsEval = newPtsEval(onPath);

    newX = interp1(t,path(:,1),newPtsEval);
    newY = interp1(t,path(:,2),newPtsEval);
    newPoints = [newX,newY];
end