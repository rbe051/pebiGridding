function [newPoints] = interpLine(path, dt)
    distS = sqrt(sum(diff(path,[],1).^2,2));
    t = [0; cumsum(distS)];
    
    newPtsEval = [0; cumsum(dt)];
    newPtsEval(end) = t(end); % Last point can not move

    newX = interp1(t,path(:,1),newPtsEval);
    newY = interp1(t,path(:,2),newPtsEval);
    newPoints = [newX,newY];
end