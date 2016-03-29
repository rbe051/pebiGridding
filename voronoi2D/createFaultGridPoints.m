function [F] = createFaultGridPoints(F,faultGridSize,circleFactor,fCut,fwCut) 
% Places fault grid points on both sides of a fault
% Arguments:
%   faultLine       k*n array of points, [x,y], describing the fault
%   fracDs          Desired distance between fault points
%   circleFactor    ratio between fracDs and circles used to create
%                   points
%
% varargin:
%   distFunc        Function setting the grid spacing
%
% Returns:
%   Pts             Fault points
%   gridSpacing     Grid spacing for each fault point
%   circCenter      Center of each circle used for creating the fault
%                   points
%   circRadius      The radius of the above circles
%   f2c             Mapping from fault points to circles.
%   c2f             Mapping from circles to fault points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:F.lines.nFault  % create fault points
  faultLine      = F.lines.lines{i};
  sePtn         = .5*[fwCut(i)==2|fwCut(i)==3; fwCut(i)==1|fwCut(i)==3];
  [p, fracSpace, fCi, fRi, f2ci,cPos, c2fi,fPos] =   ...
                          faultLinePts(faultLine,     ... 
                                       faultGridSize,...
                                       circleFactor, ...
                                       fCut(i),sePtn);
  nl = size(p,1)/2;
  if nl==0 % No fault points created
    F.lines.faultPos = [F.lines.faultPos; F.lines.faultPos(end)];
    continue
  end
  F.f.Gs          = [F.f.Gs;fracSpace];
  F.lines.faultPos = [F.lines.faultPos; size(F.f.pts,1)+1+size(p,1)];
  F.f.pts         = [F.f.pts;p];
  F.c.CC     = [F.c.CC; fCi];
  F.c.R      = [F.c.R; fRi]; 
  F.f.c         = [F.f.c; f2ci + size(F.c.fPos,1)-1];
  F.f.cPos      = [F.f.cPos; cPos(2:end) + F.f.cPos(end)-1];
  F.c.f      = [F.c.f; c2fi + size(F.f.pts,1)-nl*2];
  F.c.fPos   = [F.c.fPos; fPos(2:end) + F.c.fPos(end)-1];
end

% Add well-fault crossings
% write this as a function
endCirc = F.f.c(F.f.cPos(F.lines.faultPos([false;fwCut==1|fwCut==3]))-1);
strCirc = F.f.c(F.f.cPos(F.lines.faultPos(fwCut==2|fwCut==3)));
p       = circCircInt(F.c.CC(strCirc,:), F.c.R(strCirc),...
                      F.c.CC(endCirc,:), F.c.R(endCirc));
fId = (size(F.f.pts,1)+1:size(F.f.pts,1) + size(p,1))';
fId = reshape(fId,2,[]);
fId = repmat(fId,2,1);
cId = reshape([strCirc,endCirc]',[],1);
c2fId = repmat(F.c.fPos(cId)',2,1);
c2fId = c2fId(:);

F.f.pts = [F.f.pts;p];
nGs   = repmat(sqrt(sum(diff(p).^2,2)),1,2)';
F.f.Gs  = [F.f.Gs;reshape(nGs(:,1:2:end),[],1)];
F.c.fPos = F.c.fPos + ...
  cumsum(accumarray([cId+1;size(F.c.fPos,1)],2*[ones(1,size(cId,1)),0]));
F.c.f = insertVec(F.c.f, fId(:), c2fId);
cId = repmat(reshape(cId,2,[]),2,1);
F.f.cPos = [F.f.cPos;F.f.cPos(end)+2*cumsum(ones(size(p,1),1))];
F.f.c = [F.f.c;cId(:)];



% Remove duplicate fault Centers
if ~isempty(F.f.pts)
  [~, IA, IC] = uniquetol(F.c.CC,'byRows',true);
  F.c.CC      = F.c.CC(IA,:);
  F.c.R       = F.c.R(IA);
  [~,I]       = sort(IC);
  map         = [F.c.fPos(1:end-1), F.c.fPos(2:end)-1];
  map         = map(I,:);
  map         = arrayfun(@colon, map(:,1),map(:,2),'uniformOutput',false);
  F.c.f = F.c.f(cell2mat(map'));
  fNum        = diff(F.c.fPos);
  F.c.fPos      = cumsum([1; accumarray(IC,fNum)]);
  F.f.c = IC(F.f.c);
  % Merge intersections
  [F] = fixIntersections(F);
%                      [faultPts, fR, f2c, f2cPos, center2fault,c2fPos] =...
%       fixIntersections(faultPts, fC, fR, ...
%                        f2c, f2cPos, c2f, c2fPos);
end
   
end

function [Pts, gridSpacing, circCenter, circRadius, f2c,f2cPos, c2f,c2fPos] = ...
    faultLinePts(faultLine, fracDs, circleFactor, isCut,sePtn, varargin) 

    % Load options
    fh  = @(x) fracDs*constFunc(x);
    opt = struct('distFunc', fh);
    opt = merge_options(opt, varargin{:});
    fh = opt.distFunc;
    assert(0.5<circleFactor && circleFactor<1)
    assert(size(faultLine,1)>1 && size(faultLine,2)==2);
    
    % interpolate fault line to get desired grid spacing. 
    circCenter = interFaultLine(faultLine, fh, fracDs,sePtn);
    numOfFracPts = size(circCenter,1)-1;
    
    % Test faultLine
    if numOfFracPts == 1
      d = sqrt(sum((circCenter(2,:)-circCenter(1,:)).^2, 2));
      if d < 0.8*fh((circCenter(2,:)+circCenter(1,:))/2);
        Pts         = [];
        gridSpacing = [];
        circCenter  = [];
        circRadius  = [];
        f2c         = [];
        f2cPos      = [];
        c2f         = [];
        c2fPos      = [];
        return
      end
    end
    % Calculate the line lenth and circle radiuses. If you experience
    % imaginary faultOffset you might want to try the max lineLength
    % instead of the mean.
    lineLength = sqrt(sum((circCenter(2:end,:)-circCenter(1:end-1,:)).^2, 2));
    circRadius = circleFactor*[lineLength(1); ...
                              (lineLength(1:end-1) + lineLength(2:end))/2; ...
                               lineLength(end)];
                             
    switch isCut
      case 1
        circRadius(end) = fh(circCenter(end,:))*circleFactor;
      case 2
        circRadius(1)   = fh(circCenter(1,:))*circleFactor;
      case 3
        circRadius(1)   = fh(circCenter(1,:))*circleFactor;
        circRadius(end) = fh(circCenter(end,:))*circleFactor;
    end
    
    % Calculate the crossing of the circles
    bisectPnt = (lineLength.^2 - circRadius(2:end).^2 + circRadius(1:end-1).^2)...
                ./(2*lineLength);
    faultOffset = sqrt(circRadius(1:end-1).^2 - bisectPnt.^2);
    n1 = (circCenter(2:end,:)-circCenter(1:end-1,:))./repmat(lineLength,1,2); %Unit vector
    n2 = [-n1(:, 2), n1(:,1)];                                                %Unit normal
    
    % Set fault points on left and right side of fault
    left   = circCenter(1:end-1,:) + bsxfun(@times, bisectPnt, n1)  ...
             + bsxfun(@times, faultOffset, n2);
    right  = circCenter(1:end-1,:) + bsxfun(@times, bisectPnt, n1)  ...
             - bsxfun(@times, faultOffset, n2);
         
    % Put together result
    Pts = [right;left];
    f2c = cumsum(accumarray((1:2:size(Pts,1)+1)',1));
    f2c = repmat(f2c(2:end),2,1); 
    f2cPos = (1:2:numel(f2c)+1)';
    nf  = size(left,1);
    c2f = [   nan,          nan,         1,       nf+1;...
           (1:nf-1)', (nf+1:2*nf-1)', (2:nf)', (nf+2:2*nf)';...
              nf,           2*nf,       nan,     nan]';
    c2f = c2f(3:end-2)';
    c2fPos = [1;(3:4:numel(c2f))';numel(c2f)+1];
    gridSpacing = 2*[faultOffset;faultOffset];
end

function [p] = interFaultLine(line, fh, lineDist,sePtn, varargin)
    % Interpolate a fault line. 
    % Arguments:
    %   line        Coordinates of the fault line. Must be ordered.
    %   fh          A function handle for the relative distance function 
    %               for the interpolation fh = 1 will give a equiv distant
    %               interpolation
    %   lineDist    Scalar which set the distance between interpolation
    %               points (Relative to fh = 1)
    % varargin:    
    %               Arguments passed to fh

    % Parameters
    TOL = 1e-4; maxIt = 10000;

    % Create initial points, equally distributed.
    p = eqInterpret(line, lineDist,sePtn);
    % add auxillary points
    if sePtn(1)~=0, p = [line(1,:);p]; end
    if sePtn(2)~=0, p = [p;line(end,:)]; end
    count=0;
    while count<maxIt
      count = count+1;
      % Calculate distances, and wanted distances
      d = distAlLine(line, p);
      pmid = (p(1:end-1,:) + p(2:end,:))/2;
      dw = fh(pmid,varargin{:});
      if sePtn(1)~=0, dw(1) = dw(1).*sePtn(1);end
      if sePtn(2)~=0, dw(end) = dw(end).*sePtn(2);end

      % Possible insert or remove points
      if sum(d - dw)>min(dw)
          [~, id] = max(d - dw);
          p = [p(1:id,:); pmid(id,:); p(id+1:end,:)];
          continue
      elseif sum(d - dw)<-max(dw)
          [~, id] = min(d - dw);
          if id == 1, id = 2; end
          p = p([1:id-1,id+1:end],:);
          continue
      end
      % If we only have external nodes, we can do nothing.
      if size(p,1)<=2, return, end
      % Move points based on desired length
      Fb = dw - d;                       % Bar forces
      Fn = Fb(1:end-1) - Fb(2:end);      % Force on internal nodes
      moveNode = Fn*0.2;                 % Movement of each internal node.
      d = d + [moveNode(1); moveNode(2:end) - moveNode(1:end-1); -moveNode(end)];
      p = interpLine(line,d);            % Update node positions

      % Terminate if Nodes have moved (relative) less  than TOL
      if all(abs(moveNode)<TOL*lineDist), break; end
    end
    
    if sePtn(1)~=0, p = p(2:end,:);end
    if sePtn(2)~=0, p = p(1:end-1,:);end
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
        distA = eucDist(lineStart, p) + eucDist(p,lineEnd);
        distB = eucDist(lineStart,lineEnd);
        indx  = find(abs(distA - distB) < TOL*distB);
        if numel(indx)==0 
            jointDist = jointDist + eucDist(line(i,:), line(i+1,:));
            continue
        elseif numel(indx)>=2
            d(indx(1:end-1)) = sqrt(sum((p(indx(1:end-1),:) ... 
                             - p(indx(2:end),:)).^2,2));
            
        end
        if indx(1)>1 && eucDist(line(i,:),p(indx(1),:))>TOL
            d(indx(1)-1) = jointDist + eucDist(line(i,:), p(indx(1),:));
        end
        jointDist = eucDist(p(indx(end),:), line(i+1,:));
    end
end


function [d] = eucDist(a, b)
    d = sqrt(sum((a - b).^2,2));
end


function [newPoints] = interpLine(path, dt)
    distS = sqrt(sum(diff(path,[],1).^2,2));
    t = [0; cumsum(distS)];
    
    newPtsEval = [0; cumsum(dt)];
    newPtsEval(end) = t(end); % Last point can not move

    newX = interp1(t,path(:,1),newPtsEval);
    newY = interp1(t,path(:,2),newPtsEval);
    newPoints = [newX,newY];
end



function [F] = fixIntersections(F)
  TOL = 10*eps;
  assert(all(diff(F.f.cPos)==2),'all points must be created from exactly 2 circles');
  
  % Find conflict circles
  I = conflictCircles(F.f.pts, F.c.CC, F.c.R);
  
  circ    = find(~cellfun(@isempty,I));
  if isempty(circ)
    return
  end
  circNum = cellfun(@numel, I(circ));
  id      = zeros(sum(circNum),1);
  circPos = cumsum([1; circNum]);
  id(circPos(1:end-1)) = 1;
  id      = cumsum(id);
  circ    = [circ(id), F.f.c([F.f.cPos(vertcat(I{circ})),...
             F.f.cPos(vertcat(I{circ}))+1])];
  
  % Find shared circle
  [neigh,neighPos] = findNeighbors(circ(:,1),F.c.f,F.c.fPos, F.f.c,F.f.cPos);
  assert(all(diff(neighPos)==2));
  
  neigh  = reshape(neigh,2,[])';
  shared = 2*any(bsxfun(@eq, neigh, circ(:,2)),2) ...
          +3*any(bsxfun(@eq, neigh, circ(:,3)),2);
  keep   = find(shared);
  circ   = circ(keep,:);
  swap   = shared(keep)==2;                % Set shared circle at third row
  circ(swap,:) = [circ(swap,1),circ(swap,3),circ(swap,2)];
  
  % Remove duplicate pairs
  [~,IA] = unique(sort(circ,2),'rows');
  circ   = circ(IA,:);

  % Calculate new radiuses
  line = [F.c.CC(circ(:,3),:),reshape(mean(reshape(F.c.CC(circ(:,1:2)',:),2,[]),1),[],2)];
  int  = lineCircInt(F.c.CC(circ(:,3),:),F.c.R(circ(:,3)), line);

  % set radius to smallest
  R = sqrt(sum((F.c.CC(circ(:,1:2),:)-[int;int]).^2,2));
  I = false(size(circ(:,1:2)));
  for i = 1:numel(R)
    if R(i)<CR(circ(i))
      F.c.R(circ(i)) = R(i);
      I(circ(:,1:2)==circ(i)) = false;
      I(i) = true;
    elseif R(i)==F.c.R(circ(i))
      I(i) = true;
    end
  end
  c = unique(circ(:,1:2));

  % Calculate new Pts
  map = arrayfun(@colon,F.c.fPos(c),F.c.fPos(c+1)-1,'uniformOutput',false)';
  fId = F.c.f(horzcat(map{:})');
%   fId = unique(fId(:));
%   fId = fId(~isnan(fId));

  [neigh,neighPos] = findNeighbors(c, F.c.f,F.c.fPos, F.f.c,F.f.cPos); % I do this twice. hmm
  assert(all(diff(neighPos)==2));
  neigh = reshape(neigh,2,[])';
  
  p = circCircInt(F.c.CC(c,:), F.c.R(c),...
                 reshape(F.c.CC(neigh',:)',4,[])',reshape(F.c.R(neigh),[],2));
  assert(isreal(p),'Failed to merge fault crossings. Possible too sharp intersections');
  F.f.pts(fId,:) = p;
  map = [F.f.cPos(fId),F.f.cPos(fId)+1]';
  F.f.c(map(:)) = [c';neigh(:,1)';c';neigh(:,1)';c';neigh(:,2)';c';neigh(:,2)'];%reshape(repmat([c',c';neigh(:,1)',neigh(:,2)'],2,1),2,[]);

  [F.f.pts, ~, IC] = uniquetol(F.f.pts,'byRows',true);
  [~,I] = sort(IC);
  map = [F.f.cPos(1:end-1), F.f.cPos(2:end)-1];
  map = map(I,:);
  map = arrayfun(@colon, map(:,1),map(:,2),'uniformOutput',false);
  F.f.c = F.f.c(cell2mat(map'));
  cNum = diff(F.f.cPos);
  F.f.cPos = cumsum([1;accumarray(IC,cNum)]);
  F.c.f = IC(F.c.f);
  for i = 1:numel(c)
    f = F.c.f(F.c.fPos(c(i)):F.c.fPos(c(i)+1)-1,:);
    b = plot(F.f.pts(f,1), F.f.pts(f,2),'.','markersize',20);
    a = plot(F.c.CC(c(i),1), F.c.CC(c(i),2),'.','markersize',20);
    delete(b)
    delete(a)
  end
end

function [neigh,neighPos] = findNeighbors(c, c2f,c2fPos, f2c,f2cPos)
map   = arrayfun(@colon, c2fPos(c),c2fPos(c+1)-1,'uniformoutput',false);
pId   = cellfun(@(c) c2f(c), map,'uniformOutput',false);
neighMap = cellfun(@(c) cell2mat(arrayfun(@colon, f2cPos(c),f2cPos(c+1)-1,'uniformOutput',false)')...
                    ,pId,'uniformOutput',false);
neigh = cellfun(@(c) f2c(c),neighMap,'uniformOutput',false);
neigh = cellfun(@unique, neigh,'uniformOutput',false);
neigh = arrayfun(@(i) neigh{i}(neigh{i}~=c(i)),1:numel(neigh),'uniformOutput',false)';
neighPos = cumsum([1;cellfun(@numel, neigh)]);
neigh = vertcat(neigh{:});
end


function [p] = circCircInt(CC1, CR1, CC2,CR2)

% Expand matrices for computation
CC1 = repmat(CC1, 1,size(CC2,2)/size(CC1,2));
CR1 = repmat(CR1, 1,size(CR2,2)/size(CR1,2));
CC1 = reshape(CC1',2,[])';
CC2 = reshape(CC2',2,[])';
CR1 = reshape(CR1',1,[])';
CR2 = reshape(CR2',1,[])';

d = sqrt(sum((CC1 - CC2).^2,2));              % Distance between centers
bisectPnt = (d.^2 - CR2.^2 + CR1.^2)./(2*d);  % Mid-Point
faultOffset = sqrt(CR1.^2 - bisectPnt.^2);    % Pythagoras
n1 = (CC2-CC1)./repmat(d,1,2);                % Unit vector
n2 = [-n1(:, 2), n1(:,1)];                    % Unit normal

% Set right left and right intersection points
left   = CC1 + bsxfun(@times, bisectPnt, n1)  ...
         + bsxfun(@times, faultOffset, n2);
right  = CC1 + bsxfun(@times, bisectPnt, n1)  ...
         - bsxfun(@times, faultOffset, n2);

% Put result together
p = reshape([right,left]',2,[])';

end



function [p] = lineCircInt(CC, CR, line)
vec = line(:,3:4) - line(:,1:2);
c2l = line(:,1:2) - CC;
a   = dot(vec,vec,2);
b   = 2*dot(c2l,vec,2);
c   = dot(c2l,c2l,2) - dot(CR,CR,2);

dist    = (b.*b - 4*a.*c);
lineHit = dist>=0;
distSqr = sqrt(dist(lineHit));

%t(:,1) = -b - distSqr./(2*a); % This is the intersection on wrong side.
t = -b + distSqr./(2*a);
p = bsxfun(@times,vec,t) + line(:,1:2);


end


