clear; close all
% addpath(../../../distmesh)
%% set gridding parameters
% domain:
% set domain function
x = [1,1];
rectangle = [0,0; x(1),x(1)];   
fd = @(p) drectangle(p, 0, x(1), 0, x(2));
% Set fixed points
corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];

% Set density function
hd = @(p) 1.2 +5*(p(:,1).^2 + p(:,2).^2);
% set grid size
gs = 0.05;

%% Create grid
[Pts,t,sorting] = distmesh2d(fd, hd, gs, rectangle, corners);

G = triangleGrid(Pts,t);
Gp = pebi(G);
%% Plot Grid
figure(2)
plotGrid(G,'facecolor','none')
axis equal tight off
figure(3)
plotGrid(Gp,'facecolor','none')
axis equal tight off
%% Save grid
figure(1)
print('../../../../master/thesis/fig/ch03/optimalDelaunayInit','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch03/optimalDelaunay','-depsc')
figure(3)
print('../../../../master/thesis/fig/ch03/optimalDelaunayPebi','-depsc')