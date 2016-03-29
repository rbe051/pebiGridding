function G = compositePebiGrid(resGridSize, pdims, varargin)
% Construct a 2D composite Pebi grid. A cartesian background grid is 
% refined around faults and wells.
%
% SYNOPSIS:
%   G = compositePebiGrid(resGridSize, pdims)
%   G = compositePebiGrid(...,'Name1',Value1,'Name2',Value2,...)
%
% PARAMETERS
%   resGridSize       - Size of the reservoir grid cells, in units of
%                       meters. 
%   pdims             - Vector, length 2, [xmax, ymax], of physical size in
%                       units of meters of the computational domain. 
%
%   wellLines         - OPTIONAL.
%                       Default value empty. A struct of vectors. Each 
%                       vector, size nw x 2, is the coordinates of a 
%                       well-trace. The well is assumed to be linear 
%                       between the coorinates. If the vector only contains 
%                       one coordinate, the well is treated as a point well.
%   wellGridFactor    - OPTINAL.
%                       Default value is 0.5. This gives the relative grid
%                       size of the well grid cells compared to reservoir 
%                       grid cells. If wellGridFactor=0.5 the well cells 
%                       will be about half the size of the reservoir cells.
%   mlqtMaxLevel      - OPTIONAL.
%                       Default value 0. Number of refinement steps around 
%                       wells. 
%   mlqtLevelSteps    - OPTIONAL.  
%                       Default value -1. Vector of length mlqtMaxLevel 
%                       which specify the radius of each refinement level.
%                       The default value -1 calls the default level step
%                       in the mlqt function.
%   faultLines        - OPTINAL
%                       Default value empty. A struct of vectors.  Each 
%                       vector, size nf x 2, is the coordinates of a 
%                       fault-trace. The fault is assumed to be linear 
%                       between the coorinates
%   faultGridFactor   - OPTINAL.
%                       Default value is 0.5. This gives the relative grid
%                       size of the fault grid cells compared to reservoir 
%                       grid cells. If faultGridFactor=0.5 the fault cells 
%                       will be about half the size of the reservoir cells.
%   circleFactor      - OPTIONAL.
%                       Default value 0.6.  Valid values are between 0.5 
%                       and 1. circleFactor controll the size of the 
%                       circles used to create the fault grid points. The 
%                       circleFactor is the ratio between the radius, 
%                       and distace between the circles. A small value will
%                       place the fault points close the the faults, while
%                       a large value will place the far from the faults.
%   fullFaultEdge     - OPTIONAL.
%                       Default value FALSE. If TRUE any points that
%                       violate the sufficient fault condition will be
%                       removed. This will guarantee that faults is traced
%                       by edes in the  grid.
%
% RETURNS:
%   G                - Valid grid definition.  
%                        The fields
%                          G.cells.tag  - is TRUE for all well cells.
%                          G.faces.tag  - is TRUE for all fault edges. 
%
% EXAMPLE:
%   fl = {[0.2,0.2;0.8,0.8]};
%   wl = {[0.2,0.8;0.8,0.2]};
%   G  = compositePebiGrid(1/10,[1,1],'wellLines',wl,'faultLines',fl)
%   cla, plotGrid(G)

  %{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

    % Set options
    opt = struct('wellLines',       {{}}, ...
                 'wellGridFactor',  0.5,  ...
                 'mlqtMaxLevel',    0,    ...
                 'mlqtLevelSteps',  -1,   ...
                 'faultLines',      {{}}, ...
                 'faultGridFactor', 0.5,  ...
                 'circleFactor',    0.6,  ...
                 'fullFaultEdge',   0);         
    
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
    [wellLines, wfCut]  = splitFaults(wellLines, opt.faultLines);
    [faultLines, fwCut] = splitFaults(faultLines, opt.wellLines);
    F.lines.nFault      = numel(faultLines);
    nWell       = numel(wellLines);

    fCut        = [zeros(nWell,1); fCut];

    % Initialize variables.
    wellType  = logical([]);  % To keep track of well points 
    wGs       = [];           % Allowed grid spacing for each point
    F.f.Gs      = [];
    wellPts   = [];           % Points
    F.f.pts     = [];
    F.c.CC = [];         % Center of circle used to create fault pts
    F.c.R  = [];         % Radius of the circle
    F.f.c     = [];         % Map from a fault to the circle center
    F.f.cPos  = 1;
    F.c.f  = [];         % Map from the circle center to a fault  
    F.c.fPos = 1;
    F.lines.faultPos  = 1;         % Start and end index of a fault in Pts
    F.lines.lines = faultLines;

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
    wGs     = [wGs; wellSpace(keep)];
    wellPts = [wellPts;p(keep,:)];

  end
  [wellPts,IA] = uniquetol(wellPts,50*eps,'byRows',true);
  wGs = wGs(IA);
  % Create fault grid points
  F = createFaultGridPoints(F, faultGridSize, circleFactor, fCut, fwCut);

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

  % Remove Conflic Points
  resPts = removeConflictPoints2(resPts, wellPts,  wGs);
  resPts = removeConflictPoints2(resPts, F.f.pts, F.f.Gs);
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
  end


  % Create Grid
  Pts = [F.f.pts;wellPts; resPts];
  %[Pts,IA] = uniquetol(Pts,50*eps,'byRows',true);
  
  G = triangleGrid(Pts);
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
  if ~isempty(F.f.pts)
    N      = G.faces.neighbors + 1;
    f2cPos = [1;F.f.cPos; F.f.cPos(end)*ones(size(Pts,1)-size(F.f.pts,1),1)];
    map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
    map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
    G.faces.tag = cellfun(@(c1,c2) numel(intersect(F.f.c(c1),F.f.c(c2)))>1, map1,map2);
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













