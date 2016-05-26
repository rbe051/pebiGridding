clear; close all

%%
ds = 1/10;

rectangle = [0,0; 1,1];   
% Add wells an faults as fixed points
fixedPts = [0,0;0,1;1,1;1,0];
uni = @(p,varargin) 2*ones(size(p,1),1);
bdr = @(p) drectangle(p,0,1,0,1);
[Pts,~,sorting] = distmesh2d(bdr ,uni, ds, rectangle, fixedPts);

G = triangleGrid(Pts);

%% Plott
plotGrid(G,'facecolor','none')