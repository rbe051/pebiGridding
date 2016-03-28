function G = compositePebiGrid(resGridSize, pdims, varargin)
    % Creates a PEBI grid adapting to faults and well traces.
    %
    % Argumets:
    %   resGridSize         Size of the reservoir grid cells
    %   pdims               [xmax, ymax], array with the size of the square
    %                       to be gridded.
    %
    % Varargin:
    %   wellLines           A struct of arrays. Each array is the 
    %                       coordinates of a well trace. If the
    %                       array only contains one coordinate, the well is
    %                       treated as a point well.
    %   wellGridFactor      The relative grid size of the well grid cells
    %                       compared to reservoir grid cells
    %   mlqtMaxLevel        Number of refinement steps to be used towards
    %                       wells
    %   mlqtLevelSteps      Array of size mlqtMaxLevel which specify the
    %                       radius of each refinement level
    %   faultLines          A struct of arrays. Each array is the
    %                       coordinates of a fault trace. 
    %   faultGridFactor     The relative grid size of the fault grid cells
    %                       compared to the reservoir grid cells
    %   circleFactor        Set the relative radius of the circles used to
    %                       create the fault cells.
    %   priOrder            Array of length number of wells + number of 
    %                       faults. Sets the priority of well and fault 
    %                       traces. First element set the priority of the
    %                       first well, last element set the priority of
    %                       the last fault.
    %   fullFaultEdge       Set to true if you wish to guarantee the faults
    %                       to be traced by edges in the PEBI grid
    %
    % Returns:
    %   G                   A mrst grid structure. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % Set options
    opt = struct('wellLines',       {{}}, ...
                 'wellGridFactor',  0.5,  ...
                 'mlqtMaxLevel',    0,    ...
                 'mlqtLevelSteps',  -1,   ...
                 'faultLines',      {{}}, ...
                 'faultGridFactor', 0.5,  ...
                 'circleFactor',    0.6,  ...
                 'fullFaultEdge',   0,    ...
                 'priOrder',        []);         
    
    opt = merge_options(opt, varargin{:});
    circleFactor = opt.circleFactor;
    
    % Set grid sizes
    wellGridSize   = resGridSize*opt.wellGridFactor;
    faultGridSize  = resGridSize*opt.faultGridFactor;
    mlqtMaxLevel   = opt.mlqtMaxLevel;
    mlqtLevelSteps = opt.mlqtLevelSteps;
    
    % Test input
    assert(resGridSize>0);
    assert(numel(pdims)==2);
    assert(all(pdims>0 ));
    assert(wellGridSize>0);
    assert(mlqtMaxLevel>=0);
    assert(faultGridSize>0);
    assert(0.5<circleFactor && circleFactor<1);
    
    % Load faults and Wells
    faultLines          = opt.faultLines;
    wellLines           = opt.wellLines;
    [faultLines, fCut]  = splitFaults(faultLines, faultLines);
    [wellLines,  ~]     = splitFaults(wellLines, wellLines);
    [wellLines, wfCut]  = splitFaults(wellLines, faultLines);
    [faultLines, fwCut] = splitFaults(faultLines, wellLines);
    nFault      = numel(faultLines);
    nWell       = numel(wellLines);
    
    isWell      = [true(nWell,1); false(nFault,1)];
    fCut        = [zeros(nWell,1); fCut];

    % Initialize variables.
    faultType   = [];           % To keep track of fault points
    wellType    = logical([]);  % To keep track of well points 
    priIndex    = [];           % Priority of each point
    gridSpacing = [];           % Allowed grid spacing for each point

    wellPts       = [];         % Points
    faultPts      = [];
    fC   = [];         % Center of circle used to create fault pts
    fR   = [];         % Radius of the circle
    f2c  = [];         % Map from a fault to the circle center
    f2cPos        = 1;
    c2f  = [];         % Map from the circle center to a fault  
    c2fPos        = 1;
    faultPos      = 1;         % Start and end index of a fault in Pts


   % Create fault and well points

  for i = 1:nWell  % create well points
    wellLine       = wellLines{i};
    [p, wellSpace] = createWellGridPoints(wellLine, wellGridSize);
    keep = 1:size(p,1);
    switch wfCut(i)
      case 1
        keep = keep(1:end-1);
      case 2
        keep = keep(2:end);
      case 3
        keep = keep(2:end-1);
    end
    gridSpacing = [gridSpacing; wellSpace(keep)];
    wellPts     = [wellPts;p(keep,:)];

  end
  [wellPts,IA] = uniquetol(wellPts,50*eps,'byRows',true);
  gridSpacing = gridSpacing(IA);

  for i = 1:nFault  % create fault points
    fracLine      = faultLines{i};
    sePtn         = .5*[fwCut(i)==2|fwCut(i)==3; fwCut(i)==1|fwCut(i)==3];
    [p, fracSpace, fCi, fRi, f2ci,cPos, c2fi,fPos] =         ...
                            createFaultGridPoints(fracLine,     ... 
                                                  faultGridSize,...
                                                  circleFactor, ...
                                                  fCut(i),sePtn);
    nl = size(p,1)/2;
    if nl==0 % No fault points created
        continue
    end
    gridSpacing = [gridSpacing;fracSpace];

    faultPos    = [faultPos; size(faultPts,1)+1+size(p,1)];
    faultPts    = [faultPts;p];
    fC          = [fC; fCi];
    fR          = [fR; fRi]; 
    f2c         = [f2c;f2ci + size(c2fPos,1)-1];
    f2cPos      = [f2cPos; cPos(2:end) + f2cPos(end)-1];
    c2f         = [c2f; c2fi + size(faultPts,1)-nl*2];
    c2fPos      = [c2fPos; fPos(2:end) + c2fPos(end)-1];
  end
  
  % Add well-fault crossings
  % write this as a function
  endCirc = f2c(f2cPos(faultPos([false;fwCut==1|fwCut==3]))-1);
  strCirc = f2c(f2cPos(faultPos(fwCut==2|fwCut==3)));
  p =circCircInt(fC(strCirc,:), fR(strCirc),...
                 fC(endCirc,:), fR(endCirc));
  fId = (size(faultPts,1)+1:size(faultPts,1) + size(p,1))';
  fId = repmat(fId',2,1);
  fId = fId(:);
  cId = reshape([strCirc,endCirc]',[],1);
  c2fId = repmat(c2fPos(cId)',2,1);
  c2fId = c2fId(:);
  
  faultPts = [faultPts;p];
  nGs      = repmat(sqrt(sum(diff(p).^2,2)),1,2)';
  gridSpacing = [gridSpacing;reshape(nGs(:,1:2:end),[],1)];
  f2cPos   = [f2cPos;f2cPos(end)+2*cumsum(ones(size(p,1),1))];
  f2c = [f2c;cId];
  c2fPos   = c2fPos + cumsum(accumarray([cId+1;size(c2fPos,1)],2*[ones(1,size(cId,1)),0]));
  c2f = insertVec(c2f, fId, c2fId);
  
 
  
  % Remove duplicate fault Centers
  if ~isempty(faultPts)
    [fC, IA, IC] = uniquetol(fC,'byRows',true);
    fR = fR(IA);
    [~,I]       = sort(IC);
    map         = [c2fPos(1:end-1), c2fPos(2:end)-1];
    map         = map(I,:);
    map         = arrayfun(@colon, map(:,1),map(:,2),'uniformOutput',false);
    c2f = c2f(cell2mat(map'));
    fNum        = diff(c2fPos);
    c2fPos      = cumsum([1; accumarray(IC,fNum)]);
    f2c = IC(f2c);
    % Merge intersections
    [faultPts, fR, f2c, f2cPos, center2fault,c2fPos] =...
        fixIntersections(faultPts, fC, fR, ...
                         f2c, f2cPos, c2f, c2fPos);
  end
  % Create reservoir grid
  dx = pdims(1)/ceil(pdims(1)/resGridSize);
  dy = pdims(2)/ceil(pdims(2)/resGridSize);
  vx = 0:dx:pdims(1);
  vy = 0:dy:pdims(2);

  [X, Y] = meshgrid(vx, vy);

  resPtsInit = [X(:), Y(:)]; % Reservoir grid points before refinement


  % Refine reservoir grid
  if ~isempty(wellPts)
      res = {};
      varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', mlqtLevelSteps};
      for i = 1:size(resPtsInit,1)
          res = [res; mlqt(resPtsInit(i,:), wellPts, resGridSize, varArg{:})];
      end
      resPts      = vec2mat([res{:,1}],2);
      resGridSize = 0.5*[res{:,2}]';
  else
      resPts      = resPtsInit;
      resGridSize = repmat(0.5*min(dx,dy),size(resPts,1),1);
  end
  gridSpacing = [gridSpacing; resGridSize];

  %f2cPos = [f2cPos;f2cPos(end)*ones(size(resPts,1),1)];
  Pts       = [wellPts;faultPts;resPts];
  priIndex  = [zeros(size(wellPts,1),1);zeros(size(faultPts,1),1);ones(size(resPts,1),1)];
  wellType  = [true(size(wellPts,1),1);false(size(faultPts,1),1); false(size(resPts,1),1)];
  faultType = [false(size(wellPts,1),1);true(size(faultPts,1),1); false(size(resPts,1),1)];

  % Remove Conflic Points
  [~, ~, removed] = removeConflictPoints(Pts, gridSpacing, priIndex, wellType);
  keep = ~removed|faultType|wellType;
  Pts = Pts(keep,:);
  faultType = faultType(keep);
  wellType = wellType(keep);
  if opt.fullFaultEdge
%         [faultCenter, faultRadius] = removeFaultCircles(faultCenter, ...
%                                                         faultRadius,...
%                                                         faultPos,...
%                                                         faultToCenter,...
%                                                         removed);
      [Pts, removed] = enforceSufficientFaultCondition(Pts, ...
                                                       fC,...
                                                       fR);
      % If you whish to use the updated faultCenter and faultRadius later 
      % you can retrive them by replacing the ~ outputs.

      faultType = faultType(~removed);
      wellType = wellType(~removed);
  end


  % Create Grid
  %[Pts,IA] = uniquetol(Pts,50*eps,'byRows',true);
  Tri = delaunayTriangulation(Pts);
  G = triangleGrid(Pts, Tri.ConnectivityList);
  G = pebi(G);
  subplot(1,2,1)
  plotGrid(G,'facecolor','none')
  axis equal off
  subplot(1,2,2)
  plotGrid(G,'facecolor','none')
  axis equal off
  hold on
  for i = 1:numel(faultLines)
    l = faultLines{i};
    plot(l(:,1),l(:,2),'r')
  end


  % label fault faces.
  if any(faultType)
    f2c = f2c;

    N      = G.faces.neighbors + 1;
    f2cPos = [1;f2cPos];
    map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
    map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
    G.faces.tag = cellfun(@(c1,c2) numel(intersect(f2c(c1),f2c(c2)))>1, map1,map2);
    %G.faces.tag = logical(ft1 == ft2 & ft1 > 0 & ft2 >0);
  end
  %Label well cells
  G.cells.tag = logical(wellType);

end

function [Pts, removed] = enforceSufficientFaultCondition(Pts, CC, CR)
  I = conflictCircles(Pts, CC, CR);
  removed = false(size(Pts,1),1);
  removed(vertcat(I{:})) = true;
  Pts = Pts(~removed,:);
end


function [splitLines,isCut] = splitFaults(L1, L2)
TOL        = 20*eps;
splitLines = cell(0);
isCut      = [];
if numel(L2)==0
  isCut = 0;
  splitLines = L1;
  return
end
for i = 1:numel(L1)
  % Pick lines
  l1      = L1{i};
  l2      = L2;%([1:i-1,i+1:end]);
  linePos = cumsum([1,cellfun(@(c) size(c,1),l2)]);
  keep    = true(linePos(end)-1,1);
  keep(linePos) = false;
  l2      = cell2mat(l2');
  l1      = [l1(1:end-1,:),l1(2:end,  :)];
  l2      = [l2(1:end-1,:),l2(2:end,  :)];
  l2      = l2(keep(2:end-1),:);
  % Compute intersections
  out     = lineSegmentIntersect(l1,l2);
  [k,j]   = find(out.intAdjacencyMatrix');
  if isempty(k)
    splitLines = [splitLines, L1(i)];
    isCut      = [isCut; 0];
    continue
  end
  newPts  = [diag(out.intMatrixX(j,k)), diag(out.intMatrixY(j,k))];
  [~,I]  = sort(sum(bsxfun(@minus, newPts,l1(1:2)).^2,2));
  newPts  = newPts(I,:);
  newPts  = repmat(newPts',2,1);
  newPts  = reshape(newPts(:),2,[])';
  id      = repmat(j',2,1);
  id      = id(:);
  % Remove duplicates
  pts     = insertVec(L1{i},newPts,id+1);
  newLine =  mat2cell(pts,diff([0;j + 2*(1:numel(j))'-1;size(pts,1)]),2)';
  arg     = {'rows';'stable'};
  arg     = repmat(arg,size(newLine));
  [newLine] = cellfun(@unique, newLine,arg(1,:),arg(2,:),'UniformOutput',false);
  numPts  = cellfun(@(c) size(c,1),newLine);
  newLine = newLine(numPts>1);
  startInt = numPts(1)==1;
  endInt   = numPts(end)==1;
 
  % Split line
  splitLines = [splitLines,newLine];    
  isCut      = [isCut;...
    [ones(numel(newLine)-1,1);endInt]+2*[startInt;ones(numel(newLine)-1,1)]];
  
  for i = 1:numel(newLine)
    plot(newLine{i}(:,1),newLine{i}(:,2))
    hold on
  end
end

end

function C = insertVec(A,B,id)
[id,I] = sort(id);
B      = B(I,:);
C = zeros(size(A)+[size(B,1),0])+nan;
C(id + (0:numel(id)-1)',:) = B;
C(isnan(C)) = A;
end


function [I] = conflictCircles(Pts, CC, CR)
    TOL = 10*eps;

    nc = size(CC,1);
    np = size(Pts,1);
    removed = zeros(np,1);
    CRSqr = CR.^2;
    I = cell(numel(CR),1);
    for i = 1:nc
        distSqr = sum((repmat(CC(i,:),np,1)-Pts).^2,2);
        I{i} = find(distSqr<CRSqr(i)-10*eps);
    end                
end


function [Pts, CR,f2c,f2cPos,c2f,c2fPos] = ...
             fixIntersections(Pts, CC, CR, f2c,f2cPos, c2f,c2fPos)
  TOL = 10*eps;
  assert(all(diff(f2cPos)==2),'all points must be created from exactly 2 circles');
  
  % Find conflict circles
  I = conflictCircles(Pts, CC, CR);
  
  circ    = find(~cellfun(@isempty,I));
  if isempty(circ)
    return
  end
  circNum = cellfun(@numel, I(circ));
  id      = zeros(sum(circNum),1);
  circPos = cumsum([1; circNum]);
  id(circPos(1:end-1)) = 1;
  id      = cumsum(id);
  circ    = [circ(id), f2c([f2cPos(vertcat(I{circ})),...
             f2cPos(vertcat(I{circ}))+1])];
  
  % Find shared circle
  [neigh,neighPos] = findNeighbors(circ(:,1),c2f,c2fPos, f2c,f2cPos);
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
  line = [CC(circ(:,3),:),reshape(mean(reshape(CC(circ(:,1:2)',:),2,[]),1),[],2)];
  int  = lineCircInt(CC(circ(:,3),:),CR(circ(:,3)), line);

  % set radius to smallest
  R = sqrt(sum((CC(circ(:,1:2),:)-[int;int]).^2,2));
  I = false(size(circ(:,1:2)));
  for i = 1:numel(R)
    if R(i)<CR(circ(i))
      CR(circ(i)) = R(i);
      I(circ(:,1:2)==circ(i)) = false;
      I(i) = true;
    elseif R(i)==CR(circ(i))
      I(i) = true;
    end
  end
  c = unique(circ(:,1:2));

  % Calculate new Pts
  map = arrayfun(@colon,c2fPos(c),c2fPos(c+1)-1,'uniformOutput',false)';
  fId = c2f(horzcat(map{:})');
%   fId = unique(fId(:));
%   fId = fId(~isnan(fId));

  [neigh,neighPos] = findNeighbors(c, c2f,c2fPos, f2c,f2cPos); % I do this twice. hmm
  assert(all(diff(neighPos)==2));
  neigh = reshape(neigh,2,[])';
  
  p = circCircInt(CC(c,:), CR(c),...
                 reshape(CC(neigh',:)',4,[])',reshape(CR(neigh),[],2));
  assert(isreal(p),'Failed to merge fault crossings. Possible too sharp intersections');
  Pts(fId,:) = p;
  map = [f2cPos(fId),f2cPos(fId)+1]';
  f2c(map(:)) = [c';neigh(:,1)';c';neigh(:,1)';c';neigh(:,2)';c';neigh(:,2)'];%reshape(repmat([c',c';neigh(:,1)',neigh(:,2)'],2,1),2,[]);

  [Pts, ~, IC] = uniquetol(Pts,'byRows',true);
  [~,I] = sort(IC);
  map = [f2cPos(1:end-1), f2cPos(2:end)-1];
  map = map(I,:);
  map = arrayfun(@colon, map(:,1),map(:,2),'uniformOutput',false);
  f2c = f2c(cell2mat(map'));
  cNum = diff(f2cPos);
  f2cPos = cumsum([1;accumarray(IC,cNum)]);
  c2f = IC(c2f);
  for i = 1:numel(c)
    f = c2f(c2fPos(c(i)):c2fPos(c(i)+1)-1,:);
    b = plot(Pts(f,1), Pts(f,2),'.','markersize',20);
    a = plot(CC(c(i),1), CC(c(i),2),'.','markersize',20);
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












