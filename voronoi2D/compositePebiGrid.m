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
    

    %% Set options
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
    
    % Load faults and Wells
    faultLines  = opt.faultLines;
    [faultLines, isCut] = splitFaults(faultLines);
    wellLines   = opt.wellLines;
    nFault      = numel(faultLines);
    nWell       = numel(wellLines);
    linesToGrid = [wellLines, faultLines{:}];
    isWell      = [true(nWell,1); false(nFault,1)];
    isCut       = [zeros(nWell,1); isCut];
    % Set priority index
    priOrder = opt.priOrder;
    if isempty(priOrder)
        priOrder = 1:nFault+nWell;
    else
        assert(numel(unique(priOrder)) == nFault + nWell);
    end
    % Sort faults and wells in priority
    linesToGrid = linesToGrid(priOrder);
    isWell      = isWell(priOrder);
    isCut       = isCut(priOrder);
    
    %% Test input
    assert(resGridSize>0);
    assert(numel(pdims)==2);
    assert(all(pdims>0 ));
    assert(wellGridSize>0);
    assert(mlqtMaxLevel>=0);
    assert(faultGridSize>0);
    assert(0.5<circleFactor && circleFactor<1);
    

    %% Initialize variables.
    faultType   = [];           % To keep track of fault points
    wellType    = logical([]);  % To keep track of well points 
    priIndex    = [];           % Priority of each point
    gridSpacing = [];           % Allowed grid spacing for each point

    Pts           = [];         % Points
    faultCenter   = [];         % Center of circle used to create fault pts
    faultRadius   = [];         % Radius of the circle
    faultToCenter = [];         % Map from a fault to the circle center
    centerToFault = [];         % Map from the circle center to a fault  
    faultPos      = [];         % Start and end index of a fault in Pts


    %% Create fault and well points
    lastFaultType = 0;
    lastCCid = 0;
    for i = 1:nFault + nWell  % create well points
        if isWell(i)
            wellLine      = linesToGrid{i};
            [wellPts,wellSpace] = createWellGridPoints(wellLine, wellGridSize);
            
            np            = size(wellPts,1);
            faultType     = [faultType; zeros(np,1)];
            wellType      = [wellType; true(np,1)];
            priIndex      = [priIndex; (2+i)*ones(np,1)];
            gridSpacing   = [gridSpacing; wellSpace];
            Pts           = [Pts; wellPts];
            faultToCenter = [faultToCenter; zeros(size(wellPts,1),2)];
        end
    end

    for i = 1:nFault + nWell % create fault points
        if ~isWell(i)
            fracLine      = linesToGrid{i};    
            [faultPts, fracSpace, CC, CR, f2c, c2f] =                  ...
                                    createFaultGridPoints(fracLine,     ... 
                                                          faultGridSize,...
                                                          circleFactor, ...
                                                          isCut(i));
            nl = size(faultPts,1)/2;
            if nl==0 % No fault points created
                continue
            end

            newFaultType  = lastFaultType+1:lastFaultType+nl;
            lastFaultType = newFaultType(end);
            faultType     = [faultType; newFaultType';newFaultType'];
            wellType      = [wellType; false(2*nl,1)];
            priIndex      = [priIndex; (2+i)*ones(2*nl,1)];
            gridSpacing   = [gridSpacing;fracSpace];
            if isempty(faultPos)
                faultPos = size(Pts,1) + 1;
            end
            faultPos      = [faultPos; size(Pts,1)+1+size(faultPts,1)];
            Pts           = [Pts;faultPts];
            faultCenter   = [faultCenter; CC];
            faultRadius   = [faultRadius; CR]; 
            faultToCenter = [faultToCenter;f2c + size(centerToFault,1)];
            centerToFault = [centerToFault;c2f + size(Pts,1)-nl*2];
        end
    end
    % Remove duplicate fault Centers
    [faultCenter, IA, IC] = uniquetol(faultCenter,'byRows',true);
    faultRadius = faultRadius(IA);
    centerToFault = centerToFault(IA,:);
    faultToCenter = IC(faultToCenter);
    
    % Merge intersections
    [Pts, faultRadius] = fixIntersections(Pts, faultType, faultCenter, ...
                                         faultRadius, faultToCenter,  ...
                                         centerToFault);

    %% Create reservoir grid
    dx = pdims(1)/ceil(pdims(1)/resGridSize);
    dy = pdims(2)/ceil(pdims(2)/resGridSize);
    vx = 0:dx:pdims(1);
    vy = 0:dy:pdims(2);

    [X, Y] = meshgrid(vx, vy);

    resPtsInit = [X(:), Y(:)]; % Reservoir grid points before refinement

    
    % Refine reservoir grid
    if any(wellType)
        res = {};
        varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', mlqtLevelSteps};
        for i = 1:size(resPtsInit,1)
            res = [res; mlqt(resPtsInit(i,:), Pts(wellType,:), resGridSize, varArg{:})];
        end
        resPts      = vec2mat([res{:,1}],2);
        resGridSize = 0.5*[res{:,2}]';
    else
        resPts      = resPtsInit;
        resGridSize = repmat(0.5*min(dx,dy),size(resPts,1),1);
    end
    
    faultType   = [faultType; zeros(size(resPts,1),1)];
    wellType    = [wellType; false(size(resPts,1),1)];
    priIndex    = [priIndex; max(priIndex) + ones(size(resPts, 1), 1)];
    gridSpacing = [gridSpacing; resGridSize];
    
    Pts = [Pts;resPts];
    
    %% Remove Conflic Points
    [~, ~, removed] = removeConflictPoints(Pts, gridSpacing, priIndex, wellType);
      
    Pts = Pts(~removed|faultType,:);
    wellType = wellType(~removed|faultType,:);
    faultType = faultType(~removed|faultType,:);
    if false %opt.fullFaultEdge
        [faultCenter, faultRadius] = removeFaultCircles(faultCenter, ...
                                                        faultRadius,...
                                                        faultPos,...
                                                        faultToCenter,...
                                                        removed);
        [Pts, removed, ~, ~] = enforceSufficientFaultCondition(Pts, ...
                                                               faultType,...
                                                               faultCenter,...
                                                               faultRadius);
        % If you whish to use the updated faultCenter and faultRadius later 
        % you can retrive them by replacing the ~ outputs.
        
        faultType = faultType(~removed);
        wellType = wellType(~removed);
    end
    
    
    %% Create Grid
    Pts = uniquetol(Pts,50*eps,'byRows',true);
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
    N = G.faces.neighbors + 1;
    faultType = [0; faultType];
    ft1 = faultType(N(:,1));
    ft2 = faultType(N(:,2));
    G.faces.tag = logical(ft1 == ft2 & ft1 > 0 & ft2 >0);
    
    %Label well cells
    G.cells.tag = logical(wellType);

end


function [CC, CR] = removeFaultCircles(CC, CR,faultPos,faultToCenter, remFaultPts)
    assert(size(CC,1)==size(CR,1));    
    % Create a index map of the points
    faultIndex = (1:size(faultToCenter,1))';
    remCirc = zeros(size(CC,1),1);
    for i = 1:size(faultPos,1)-1 %iterate over all faults
        % Find removed fault points from current fault
        faultPts = faultIndex(faultPos(i):faultPos(i+1)-1);
        isRemovedFaults = false(size(faultToCenter,1),1);
        isRemovedFaults(faultPts) = remFaultPts(faultPts);
        % Find the first and last circle in curren fault
        Cfrom  = faultToCenter(faultPos(i));
        Cto = faultToCenter(faultPos(i+1)-1)+1;
        % Find which circles has a removed fault point
        toRemove1 = false(size(CC,1),1);
        toRemove2 = toRemove1;
        circleToRemove1 = faultToCenter(faultIndex(isRemovedFaults));
        circleToRemove1 = unique(circleToRemove1);
        % Remove circle if a fault pnt on both sides is removed
        circleToRemove2 = circleToRemove1 + 1;
        circleToRemove1 = [circleToRemove1; Cto];
        circleToRemove2 = [Cfrom; circleToRemove2];
        toRemove1(circleToRemove1) = true;
        toRemove2(circleToRemove2) = true;
        toRemove = and(toRemove1, toRemove2);
        remCirc = remCirc + toRemove;
    end
    remCirc = logical(remCirc);
    CC = CC(~remCirc,:);
    CR = CR(~remCirc);    

end



function [splitLines,isCut] = splitFaults(faultLines)
TOL = 20*eps;
splitLines = cell(0);
isCut      = [];
for i = 1:numel(faultLines)
  % Pick lines
  l1 = faultLines{i};
  l2 = faultLines([1:i-1,i+1:end]);
  linePos = cumsum([1,cellfun(@(c) size(c,1),l2)]);
  keep = true(linePos(end)-1,1);
  keep(linePos) = false;
  l2 = cell2mat(l2');
  l1 = [l1(1:end-1,:),l1(2:end,  :)];
  l2 = [l2(1:end-1,:),l2(2:end,  :)];
  l2 = l2(keep(2:end-1),:);
  % Compute intersections
  out = lineSegmentIntersect(l1,l2);
  [k,j] = find(out.intAdjacencyMatrix');
  newPts = [diag(out.intMatrixX(j,k)), diag(out.intMatrixY(j,k))];
  newPts = repmat(newPts',2,1);
  newPts = reshape(newPts(:),2,[])';
  id = repmat(j',2,1);
  id = id(:);
  % Remove duplicates
  pts = insertVec(faultLines{i},newPts,id+1);
  newLine =  mat2cell(pts,diff([0;j + 2*(1:numel(j))'-1;size(pts,1)]),2)';
  arg = {'rows';'stable'};
  arg = repmat(arg,size(newLine));
  [newLine] = cellfun(@unique, newLine,arg(1,:),arg(2,:),'UniformOutput',false);
  numPts = cellfun(@(c) size(c,1),newLine);
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
C = zeros(size(A)+[size(B,1),0])+nan;
C(id + (0:numel(id)-1)',:) = B;
C(isnan(C)) = A;
end


function [I] = conflictCircles(Pts, isFault, CC, CR)
    TOL = 10*eps;
    assert(size(isFault,1)==size(Pts,1));

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


function [Pts, CR] = fixIntersections(Pts, isFault, CC, CR, f2c, c2f)
  TOL = 10*eps;
 
  % Find conflict circles
  I = conflictCircles(Pts, isFault, CC, CR);
  
  circ = find(~cellfun(@isempty,I));
  circNum = cellfun(@numel, I(circ));
  id = zeros(sum(circNum),1);
  circPos = cumsum([1; circNum]);
  id(circPos(1:end-1)) = 1;
  id = cumsum(id);
  %rem  = cellfun(@numel, I(circ))>1; % Ignore these cases
  %circ = circ(~rem);
  
  circ = [circ(id), f2c(vertcat(I{circ}),:)]; 
  
  % Find shared circle
  neigh = findNeighbors(circ(:,1),c2f, f2c); % OBS assumes specific ordering on faults
  
  shared = 2*any(bsxfun(@eq, neigh, circ(:,2)),2) ...
          +3*any(bsxfun(@eq, neigh, circ(:,3)),2);
  keep = find(shared);
  circ = circ(keep,:);
  swap = shared(keep)==2;
  circ(swap,:) = [circ(swap,1),circ(swap,3),circ(swap,2)];
  [~,IA] = unique(sort(circ,2),'rows');
  circ = circ(IA,:);

  % Calculate new radiuses
  line = [CC(circ(:,3),:),reshape(mean(reshape(CC(circ(:,1:2)',:),2,[]),1),[],2)];
  int = lineCircInt(CC(circ(:,3),:),CR(circ(:,3)), line);
  % set radius to smallest
  R = sqrt(sum((CC(circ(:,1:2),:)-[int;int]).^2,2));
  for i = 1:numel(R)-1
    CR(circ(i)) = min(CR(circ(i)),R(i));
  end
  % Remove duplicates
  c = unique(circ(:,1:2));
      
  % Calculate new Pts
  fId = c2f(c,:);
  fId = unique(fId(:));
  fId = fId(~isnan(fId));
  
  neigh = findNeighbors(c(:), c2f, f2c); % I do this twice. hmm
  p = circCircInt(CC(c(:),:), CR(c(:)),...
                 reshape(CC(neigh',:)',4,[])',reshape(CR(neigh),[],2));
  
  Pts(fId,:) = p;
end

function [neigh] = findNeighbors(c, c2f, f2c)
%OBS does not work on end points
pId   = c2f(c,[1,3]);
neigh = reshape(f2c(pId',:)',4,[])';
neigh = reshape(neigh(bsxfun(@ne, neigh,c)),[],2);
end


function [p] = circCircInt(CC1, CR1, CC2,CR2)

CC1 = repmat(CC1, 1,size(CC2,2)/size(CC1,2));
CR1 = repmat(CR1, 1,size(CR2,2)/size(CR1,2));
CC1 = reshape(CC1',2,[])';
CC2 = reshape(CC2',2,[])';
CR1 = reshape(CR1',1,[])';
CR2 = reshape(CR2',1,[])';

d = sqrt(sum((CC1 - CC2).^2,2));
bisectPnt = (d.^2 - CR2.^2 + CR1.^2)...
                ./(2*d);
faultOffset = sqrt(CR1.^2 - bisectPnt.^2);
n1 = (CC2-CC1)./repmat(d,1,2); %Unit vector
n2 = [-n1(:, 2), n1(:,1)];     %Unit normal

% Set fault points on left and right side of fault
left   = CC1 + bsxfun(@times, bisectPnt, n1)  ...
       + bsxfun(@times, faultOffset, n2);
right  = CC1 + bsxfun(@times, bisectPnt, n1)  ...
       - bsxfun(@times, faultOffset, n2);

% Put together result
p = [right;left];


end



function [p] = lineCircInt(CC, CR, line)
p = nan(size(CC,1),2);
t = zeros(size(CC,1),1);
vec = line(:,3:4) - line(:,1:2);
c2l = line(:,1:2) - CC;
a   = dot(vec,vec,2);
b   = 2*dot(c2l,vec,2);
c = dot(c2l,c2l,2) - dot(CR,CR,2);

dist = (b.*b - 4*a.*c);
lineHit = dist>=0;
distSqr = sqrt(dist(lineHit));

%t(:,1) = -b - distSqr./(2*a); % This is intersection on wrong side.
t = -b + distSqr./(2*a);

p = bsxfun(@times,vec,t) + line(:,1:2);


end












