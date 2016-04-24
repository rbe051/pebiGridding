function G = pebiGrid(resGridSize, pdims, varargin)
% Construct a 2D Pebi grid.
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
%   wellGridFactor    - OPTIONAL.
%                       Default value is 0.5. This gives the relative grid
%                       size of the well grid cells compared to reservoir 
%                       grid cells. If wellGridFactor=0.5 the well cells 
%                       will be about half the size of the reservoir cells.
%   wellRefinement    - OPTIONAL
%                       Default value FALSE. Set to true to turn on
%                       refinement around wells.
%   epsilon           - OPTIONAL
%                       Default value 0.25/max(pdims). epsilon set the
%                       refinement transition around wells. The density
%                       function for the reservoir grid is set by
%                       rho~exp(-distance to well / epsilon).
%   faultLines        - OPTIONAL
%                       Default value empty. A struct of vectors.  Each 
%                       vector, size nf x 2, is the coordinates of a 
%                       fault-trace. The fault is assumed to be linear 
%                       between the coorinates
%   faultGridFactor   - OPTIONAL.
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
%                       Default value TRUE. If TRUE any points that
%                       violate the sufficient fault condition will be
%                       removed. This will guarantee that faults is traced
%                       by edes in the  grid.
%
% RETURNS:
%   G                - Valid grid definition.  
%                        The fields
%                          - G.cells.tag is TRUE for all well cells.
%                          - G.faces.tag is TRUE for all fault edges. 
%
% EXAMPLE:
%   fl = {[0.2,0.2;0.8,0.8]};
%   wl = {[0.2,0.8;0.8,0.2]};
%   G  = compositePebiGrid(1/10,[1,1],'wellLines',wl,'faultLines',fl)
%   cla, plotGrid(G)
%
% SEE ALSO
%   compositePebiGrid, pebi, createFaultGridPoints, createWellGridPoints.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%
% distMesh is used to create the background grid 
% (http://persson.berkeley.edu/distmesh/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

% Set options
opt = struct('wellLines',       {{}}, ...
             'wellGridFactor',  0.5, ...
             'wellRefinement',  false, ...
             'epsilon',         -1,...
             'faultLines',      {{}}, ...
             'faultGridFactor', 0.5, ...
             'circleFactor',    0.6);  

opt          = merge_options(opt, varargin{:});
circleFactor = opt.circleFactor;
wellRef      =  opt.wellRefinement;
wellEps      = opt.epsilon;

% Set grid sizes
wellGridFactor  = opt.wellGridFactor;
faultGridFactor = opt.faultGridFactor;
wellGridSize    = resGridSize*wellGridFactor;
faultGridSize   = resGridSize*faultGridFactor;
if wellEps<0
    wellEps = 0.25/max(pdims);
end

% Test input
assert(resGridSize>0);
assert(numel(pdims)==2);
assert(all(pdims>0 ));
assert(wellGridSize>0);
assert(faultGridSize>0);
assert(0.5<circleFactor && circleFactor<1);

% Load faults and Wells
faultLines                = opt.faultLines;
wellLines                 = opt.wellLines;
[faultLines, fCut, fwCut] = splitAtInt(faultLines, wellLines);
[wellLines,  ~, wfCut]    = splitAtInt(opt.wellLines, opt.faultLines);

% Create well Points
[wellPts, ~] = createWellGridPoints(wellLines, wellGridSize,'wfCut', wfCut);

% create distance functions
if wellRef && ~isempty(wellPts)
    hres   = @(x) min((ones(size(x,1),1)/wellGridFactor), ...
                  min(1.2*exp(pdist2(x,wellPts)/wellEps),[],2));
    hfault = @(x) wellGridSize*faultGridFactor*hres(x);
else
hres   = @(p) constFunc(p)/wellGridFactor;
hfault = @(p) constFunc(p)*faultGridSize;
end

% Create fault points
F = createFaultGridPoints(faultLines, faultGridSize,'circleFactor', circleFactor,...
                          'fCut', fCut,'fwCut', fwCut, 'distFun', hfault);

% Create Reservoir grid points
% set domain function
x = pdims;
rectangle = [0,0; x(1),x(1)];   
fd = @(p) drectangle(p, 0, x(1), 0, x(2));
% Set fixed points
corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];

% Add wells an faults as fixed points
fixedPts = [F.f.pts; wellPts; corners];

[Pts,~,sorting] = distmesh2d(fd, hres, wellGridSize, rectangle, fixedPts);

% Distmesh change the order of all points. We undo this sorting.
isFault = false(size(Pts,1),1); 
isFault(1:size(F.f.pts,1)) = true;
isFault = isFault(sorting);
[~,If]   = sort(sorting(isFault));
isWell  = false(size(Pts,1),1); 
isWell(size(F.f.pts,1)+(1:size(wellPts,1)))= true;
isWell  = isWell(sorting);
[~,Iw]   = sort(sorting(isWell));
isRes   = ~isFault & ~isWell;

fPts = Pts(isFault,:);
fPts = fPts(If,:);
wPts = Pts(isWell,:);
wPts = wPts(Iw,:);
Pts  = [fPts; wPts; Pts(isRes,:)];

% Create grid
t    = delaunay(Pts);
% Fix boundary
pmid = (Pts(t(:,1),:)+Pts(t(:,2),:)+Pts(t(:,3),:))/3;% Compute centroids
t    = t(fd(pmid)<-0.001*wellGridFactor,:);    % Keep interior triangles

G = triangleGrid(Pts, t);
G = pebi(G);

% label fault faces.
if ~isempty(F.f.pts)
  N      = G.faces.neighbors + 1; 
  % N == 1 is now a boundary fault, so we have to add 1 to the start of 
  % cPos.  We also add empty mapping for no fault pts.
  f2cPos = [1;F.f.cPos; F.f.cPos(end)*ones(size(Pts,1)-size(F.f.pts,1),1)];
  map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
  map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
  G.faces.tag = cellfun(@(c1,c2) numel(intersect(F.f.c(c1),F.f.c(c2)))>1, map1,map2);
else
  G.faces.tag = false(G.faces.num,1);
end

%Label well cells
if ~isempty(wellPts)
  G.cells.tag = false(G.cells.num,1);
  % Add tag to all cells generated from wellPts
  wellCells = size(F.f.pts,1)+1:size(F.f.pts,1)+size(wellPts,1);
  G.cells.tag(wellCells)= true;
  
  % Add tag to well-fault crossings
  endOfLine = fwCut==1 | fwCut==3;        % Crossing at end of fault
  strOfLine = fwCut==2 | fwCut==3;        % Crossing at start of fault
  fIde      = F.l.fPos([false;endOfLine]);
  fIds      = F.l.fPos(strOfLine);
  fToTag    = [F.l.f(mcolon(fIde - 2,fIde-1)); F.l.f(mcolon(fIds,fIds+1))];

  G.cells.tag(fToTag) = true;
else
  G.cells.tag = false(G.cells.num,1);
end
end

