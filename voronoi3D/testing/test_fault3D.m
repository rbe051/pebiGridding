clc; clear;close all


%%
dtPts = [-1,-1,-1;  ...
          1,-1,-1;  ...
          1, 1,-1;  ...
         -1, 1,-1;  ...
         -1,-1, 1;  ...
          1,-1, 1;  ...
          1, 1, 1;  ...
         -1, 1, 1];
%% two horizontal faults
n = 10;
x = linspace(-1,1,(n));
[X,Y] = meshgrid(x,x);
facePts = [X(:),Y(:)];

face.ConnectivityList = delaunay(facePts);
face.Points = [facePts,0.33*ones(size(facePts,1),1)];
face2 = face;
face2.Points(:,3) = face2.Points(:,3) - 0.67;


rho = @(p) max(diff(x))*0.75*ones(size(p,1),1);



[faultPts,CC,CR] = createFaultGridPoints3D(face, rho);
[faultPts2,CC2,CR2] = createFaultGridPoints3D(face2, rho);

faultPts = [faultPts;faultPts2];
CC = [CC;CC2];
CR = [CR;CR2];

%% One tilted fault
n = 10;
x = linspace(-1,1,(n));
[X,Y] = meshgrid(x,x);
facePts = [X(:),Y(:)];

face.ConnectivityList = delaunay(facePts);
face.Points = [facePts,zeros(size(facePts,1),1)];
face.Points(:,3) = face.Points(:,3) + 0.5*face.Points(:,1);


rho = @(p) max(diff(x))*0.8*ones(size(p,1),1);



[faultPts,CC,CR] = createFaultGridPoints3D(face, rho);


%% Generate Grid POints and grid
gridPts = rand(n^3,3)*2-1;
%x = linspace(-1,1,n);
%[X,Y,Z] = ndgrid(x,x,x);
%gridPts = [X(:),Y(:),Z(:)];
pts = [gridPts;faultPts];

[pts,removed] = faultSufCond(pts,CC,CR);

G = voronoi3D(pts, dtPts);

plotGrid(G)
