clear; close all
addpath ../../../distmesh/
addpath ../

%%
load('showInnerOuterCellsElephant.mat');

%%
rectangle = [min(bdr); max(bdr)]; 
% Add wells an faults as fixed points
fixedPts = [0.4521,-.2955; ...
            0.04004,-0.3369;...
            0.003358, -0.4056;...
            0.1806,-0.4547;...
            0.2467, -.4492;...
            0.1795,-.1387; ...
            0.4963, -.2002;...
            .4784,-.1405;...
            .3748,-.3206];
gS = 4/(size(pts,1));
uni = @(p,varargin) 2*ones(size(p,1),1);
[Pts,t,sorting] = distmesh2d(@signDistFuncPolygon ,uni, gS, rectangle, fixedPts,'bdr',bdr);

Gt = triangleGrid(Pts, t);
G = pebi(Gt);

%%
plotGrid(G)