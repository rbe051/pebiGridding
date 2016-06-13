function [G, Gn] = pebiGridForCompareFractures(resGridSize, pdims, varargin)
% Construct a 2D Pebi grid. The differance from pebiGrid(...) is that this
% version creates two grids. One with centroids following the wells, and 
% one where edges follows the well. All other sites are the same
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
%
% RETURNS:
%   G                - Valid grid definition. With centroid centerd wells
%                        The fields
%                          - G.cells.tag is TRUE for all well cells.
%                          - G.faces.tag is TRUE for all fault edges. 
%   Gn               - Valid grid definition. Edge centered wells.
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
             'protLayer',       false, ...
             'circleFactor',    0.6, ...
             'protD',           {{@(p) ones(size(p,1),1)*resGridSize/10}});

opt          = merge_options(opt, varargin{:});
circleFactor = opt.circleFactor;
wellRef      =  opt.wellRefinement;
wellEps      = opt.epsilon;

% Set grid sizes
wellGridFactor  = opt.wellGridFactor;

wellGridSize    = resGridSize*wellGridFactor;

if wellEps<0
    wellEps = 0.25/max(pdims);
end

% Test input
assert(resGridSize>0);
assert(numel(pdims)==2);
assert(all(pdims>0 ));
assert(wellGridSize>0);
assert(0.5<circleFactor && circleFactor<1);

% Load faults and Wells
faultLines                = opt.wellLines;
wellLines                 = opt.wellLines;
[faultLines, fCut, fwCut] = splitAtInt(faultLines, {});
[wellLines,  wCut, wfCut,IC] = splitAtInt(wellLines, {});

% Load protection layer
protD = opt.protD;
if numel(protD) == 1
  protD = repmat(protD,numel(wellLines),1);num2cell(protD, 2);
else
  assert(numel(protD) == numel(opt.wellLines));
  protD = protD(IC);
end

%SE = repmat(0.5,numel(wellLines),2);
[wellPts, wGs,protPts,pGs] = createWellGridPoints(wellLines, wellGridSize,...
                                                 'wCut',wCut,'protLayer',opt.protLayer,...
                                                 'protD',protD);


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
F = createFaultGridPoints(faultLines, wellGridSize,'circleFactor', circleFactor,...
                          'fCut', fCut);

% Create Reservoir grid points
% set domain function
x = pdims;
rectangle = [0,0; x(1),x(1)];   
fd = @(p) drectangle(p, 0, x(1), 0, x(2));
% Set fixed points
corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];

% Add wells an faults as fixed points
fixedPts = [F.f.pts; corners];

[Pts,~,sorting] = distmesh2d(fd, hres, wellGridSize, rectangle, fixedPts);

% Distmesh change the order of all points. We undo this sorting.
isFault = false(size(Pts,1),1); isFault(1:size(F.f.pts,1)) = true;
isFault = isFault(sorting);
[~,If]   = sort(sorting(isFault));

isRes   = ~isFault;

fPts = Pts(isFault,:);
fPts = fPts(If,:);

Ptsn  = [fPts; Pts(isRes,:)];
Pts  = [protPts;wellPts; Pts(isRes,:)];

% Create grid With centroids following well: 
t    = delaunay(Pts);
% Fix boundary
pmid = (Pts(t(:,1),:)+Pts(t(:,2),:)+Pts(t(:,3),:))/3;% Compute centroids
t    = t(fd(pmid)<-0.001*wellGridFactor,:);    % Keep interior triangles

bdr = [0,0;0,pdims(2);pdims(1),pdims(2);pdims(2),0];
%Gt = triangleGrid(Pts,t);
%G = pebi(Gt);
G = clippedPebi2D(Pts,bdr);

%Label well cells
if ~isempty(wellPts)
  G.cells.tag = false(G.cells.num,1);
  % Add tag to all cells generated from wellPts
  wellCells = size(protPts,1)+1:size(protPts,1)+size(wellPts,1);
  G.cells.tag(wellCells)= true;
else
  G.cells.tag = false(G.cells.num,1);
end

% tag protection sites
G.cells.protectionCells = false(G.cells.num,1);
G.cells.protectionCells(1:size(protPts,1)) = true;

% Create edge centered well


% Create grid
tn    = delaunay(Ptsn);
% Fix boundary
pmidn = (Ptsn(tn(:,1),:)+Ptsn(tn(:,2),:)+Ptsn(tn(:,3),:))/3;% Compute centroids
tn    = tn(fd(pmidn)<-0.001*wellGridFactor,:);    % Keep interior triangles

%Gtn = triangleGrid(Ptsn,tn);
%Gn = pebi(Gtn);
Gn = clippedPebi2D(Ptsn,bdr);

% label fault faces.
if ~isempty(F.f.pts)
  N      = Gn.faces.neighbors + 1; 
  % N == 1 is now a boundary fault, so we have to add 1 to the start of 
  % cPos.  We also add empty mapping for no fault pts.
  f2cPos = [1;F.f.cPos; F.f.cPos(end)*ones(size(Ptsn,1)-size(F.f.pts,1),1)];
  map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
  map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
  Gn.faces.tag = cellfun(@(c1,c2) numel(intersect(F.f.c(c1),F.f.c(c2)))>1, map1,map2);
else
  Gn.faces.tag = false(Gn.faces.num,1);
end


end