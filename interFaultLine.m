function [p] = interFaultLine(line, fh, lineDist, varargin)
    % Interpolate a fault line. 
    % Arguments:
    %   line        Coordinates of the fault line. Must be ordered.
    %   fh          A function handle for the relative distance function 
    %               for the interpolation fh = 1 will give a equiv distant
    %               interpolation
    %   lineDist    Scaler which set the distance between interpolation
    %               points (Relative to fh = 1)
    %   varargin    Arguments passed to fh

    % Parameters
    TOL = 1e-4; maxIt = 10000;

    % Create initial points, equally distributed.
    p = eqInterpret(line, lineDist);

    count=0;
    while count<maxIt
      count = count+1;
      % Calculate distances, and wanted distances
      d = distAlLine(line, p);
      pmid = (p(1:end-1,:) + p(2:end,:))/2;
      dw = lineDist*fh(pmid,varargin{:}); % Multiply by lineDist since fh is 
                                          % the relative size fuction

      % Possible insert or remove points
      if sum(d - dw) > min(dw)
          id = find(min(d));
          p = [p(1:id,:); pmid(id,:); p(id+1:end,:)];
          continue
      elseif sum(d - dw) < - max(dw)
          id = find(min(d));
          if id == 1, id = 2; end
          p = p([1:id-1,id+1:end],:);
          continue
      end
      % If we only have external nodes, we can do nothing.
      if size(p,1)<=2, return, end
      % Move points based on desired length
      Fb = dw - d;                       % Bar forces (scalars)
      Fn = Fb(1:end-1) - Fb(2:end);      % Force on internal nodes
      moveNode = Fn*0.2;                 % Movement of each internal node.
      d = d + [moveNode(1); moveNode(2:end) - moveNode(1:end-1); -moveNode(end)];
      p = interpLine(line,d);            % Update node positions
        
      % Terminate if Nodes have moved (relative) less  than TOL
      if all(moveNode<TOL*lineDist), return; end
    end

    if count == maxIt
        warning('Fault interpolation did not converge.')
    end

end


function [d] = distAlLine(line, p)
% Calculates the distace between consecutive interpolation points along
% line
% Arguments:
%   line    line that is interpolated
%   p       Interpolation points
% Returns:
%   d       distance between consecutive points of p, along line
    TOL = 50*eps;
    N = size(p,1);
    d = zeros(N-1,1);
    jointDist = 0;
    for i = 1:size(line,1)-1
        lineStart = repmat(line(i,:), size(p,1),1);
        lineEnd = repmat(line(i+1,:), size(p,1),1);
        distA = eucDistSqrt(lineStart, p) + eucDistSqrt(p,lineEnd);
        distB = eucDistSqrt(lineStart,lineEnd);
        indx  = find(abs(distA - distB) < TOL); %Find points on line segment
        if numel(indx)==0 
            jointDist = jointDist + eucDistSqrt(line(i,:), line(i+1,:));
            continue
        elseif numel(indx)>=2
            d(indx(1:end-1)) = sqrt(sum((p(indx(1:end-1),:) ... 
                             - p(indx(2:end),:)).^2,2));
            
        end
        if indx(1)>1 && eucDistSqrt(line(i,:),p(indx(1),:))>TOL
            d(indx(1)-1) = jointDist + eucDistSqrt(line(i,:), p(indx(1),:));
        end
        jointDist = eucDistSqrt(p(indx(end),:), line(i+1,:));
    end
end


function [d] = eucDistSqrt(a, b)
    d = sqrt(sum((a - b).^2,2));
end





