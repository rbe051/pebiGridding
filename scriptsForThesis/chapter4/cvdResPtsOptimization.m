clear; close all
%addpath ../../voronoi2D/
%addpath ../
%addpath ../../voronoi3D/
%addpath ../../../vem/mat/
%addpath ../../../vem/mat/VEM2D
%% Create fault
f = {[0.2,0.4;0.8,0.6]};

fGs = 0.1;
F = createFaultGridPoints(f,fGs);

%% set gridding parameters
% domain:
% set domain function
x = [1,1];
% Set fixed points
bnd = [-1,-1;-1,x(2);x(1),x(2);x(1),-1];

% Set density function
rho = @(p) exp(-1*p(:,1).^2 - 1*p(:,2).^2);
% set grid size
gs = 0.02;

%% Create initial reservoir sites
n = round(norm(x)/gs);
rSi = rand(n,2);

G = createCVD(rSi,bnd,'fixedPts',F.f.pts,'rho',rho);


%% plot
plotGrid(G);
axis equal
plotLinePath(f)