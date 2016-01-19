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
      d = sqrt(sum((p(1:end-1,:) - p(2:end,:)).^2,2));
      pmid = (p(1:end-1,:) + p(2:end,:))/2;
      dw = lineDist*fh(pmid,varargin{:}); % Multiply by lineDist since fh is 
                                          % the relative size fuction

      % Possible insert or remove points
      if sum(d - dw) > min(d)
          id = find(min(d));
          p = [p(1:id,:); pmid(id,:); p(id+1:end,:)];
          continue
      elseif sum(d - dw) < - max(d)
          id = find(min(d));
          if id == 1, id = 2; end
          p = p([1:id-1,id+1:end],:);
          continue
      end

      % 3. Move points based on bar desired length
      Fb = dw - d;                       % Bar forces (scalars)
      Fn = Fb(1:end-1) - Fb(2:end);      % Force on internal nodes
      moveNode = Fn*0.2;                 % Movement of each internal node.
      d = d + [moveNode(1); moveNode(2:end) - moveNode(1:end-1); -moveNode(end)];
      p = interpLine(line,d);                           % Update node positions
        
      % Terminate if Nodes have moved (relative) less  than TOL
      if all(moveNode<TOL*lineDist), return; end
    end

    if count == maxIt
        warning('Fault interpolation did not converge.')
    end

end