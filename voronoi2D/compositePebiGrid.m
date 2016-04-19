function G = compositePebiGrid(resGridSize, pdims, varargin)
% Construct a 2D composite Pebi grid. A cartesian background grid that is 
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
             'fullFaultEdge',   1);         

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
faultLines                = opt.faultLines;
wellLines                 = opt.wellLines;
[faultLines, fCut, fwCut] = splitAtInt(faultLines, wellLines);
[wellLines,  ~, wfCut]    = splitAtInt(opt.wellLines, opt.faultLines);


% Create well points
[wellPts, wGs] = createWellGridPoints(wellLines, wellGridSize,'wfCut',wfCut);

% Create fault points
F = createFaultGridPoints(faultLines, faultGridSize, 'circleFactor', circleFactor,...
                          'fCut',fCut,'fwCut', fwCut);

% Create reservoir grid
dx = pdims(1)/ceil(pdims(1)/resGridSize);
dy = pdims(2)/ceil(pdims(2)/resGridSize);
vx = 0:dx:pdims(1);
vy = 0:dy:pdims(2);

[X, Y] = meshgrid(vx, vy);

resPtsInit = [X(:), Y(:)];

% Refine reservoir grid
if ~isempty(wellPts)
    res = {};
    varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', mlqtLevelSteps};
    for i = 1:size(resPtsInit,1)
        res = [res; mlqt(resPtsInit(i,:), wellPts, resGridSize, varArg{:})];
    end
    resPts = vec2mat([res{:,1}],2);
    %resGridSize = 0.5*[res{:,2}]';
else
    resPts = resPtsInit;
    % resGridSize = repmat(0.5*min(dx,dy),size(resPts,1),1);
end

% Remove Conflic Points
resPts = removeConflictPoints2(resPts, wellPts,  wGs);
resPts = removeConflictPoints2(resPts, F.f.pts, F.f.Gs);
resPts = removeConflictPoints2(resPts, F.c.CC, F.c.R);

% Create Grid
Pts = [F.f.pts;wellPts; resPts];
G = triangleGrid(Pts);
G = pebi(G);

% label fault faces.
if ~isempty(F.f.pts)
  N      = G.faces.neighbors + 1;
  f2cPos = [1;F.f.cPos; F.f.cPos(end)*ones(size(Pts,1)-size(F.f.pts,1),1)];
  map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
  map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
  G.faces.tag = cellfun(@(c1,c2) numel(intersect(F.f.c(c1),F.f.c(c2)))>1, map1,map2);
else
  G.faces.tag = false(G.faces.num);
end

%Label well cells
G.cells.tag = false(G.cells.num,1);
G.cells.tag(size(F.f.pts,1)+1:size(F.f.pts,1)+size(wellPts,1))= true;
G.cells.tag(F.l.fPos(end):size(F.f.pts,1)+1) = true;
end


