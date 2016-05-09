clear; close all
addpath ../
%addpath ../../../vem/mat/VEM2D
%addpath ../../../vem/mat/

%% Set boundary and parameters
k = 100;
n =     [-1, 0; ...
          0,  1; ...
          1, 0; ...
          0, -1];
x0 = 10*[-1, -1; ...
         -1,  1; ...
          1,  1; ...
          1, -1];
    

%% Create none CVD
pts = [ 0,.2; ...
       .3,.1
       .1,.8; ...
       .3,.4; ...
       .5,.6; ...
       .8,.9; ...
       .6,.3]; 

Gt = triangleGrid(pts);

G = clippedPebi2D(pts, n, x0);
G = computeGeometry(G);

%% Plotting
col = get(gca,'ColorOrder');
%% Only voronoi
figure(); hold on
plotGrid(G,'facecolor','none')
plot(pts(:,1), pts(:,2),'.','color',col(1,:),'markersize',15);
axis off equal
axis([-0.1,1,0,1])
set(gca,'position',[0 0 1 1],'units','normalized')
print('../../../../master/thesis/fig/ch02/delaunayVoronoiDualityVor','-depsc')

%% Only delaunay
figure(); hold on
plotGrid(Gt,'facecolor','none')
plot(pts(:,1), pts(:,2),'.','color',col(1,:),'markersize',15);
axis off equal
axis([-0.1,1,0,1])
set(gca,'position',[0 0 1 1],'units','normalized')
print('../../../../master/thesis/fig/ch02/delaunayVoronoiDualityDel','-depsc')


%% Both delaunay and voronoi
figure(); hold on
plotGrid(Gt,'facecolor','none')
plotGrid(G,'facecolor','none')
plot(pts(:,1), pts(:,2),'.','color',col(1,:),'markersize',20);
% Plot circumcircle
t = 5;
dt = delaunayTriangulation(pts);
[CC, r] = circumcenter(dt,t);
theta = linspace(0,2*pi)';
for i = 1:size(CC)
X = repmat(CC(i,:),100,1) + repmat(r(i),100,2).*[cos(theta), sin(theta)];
plot(X(:,1), X(:,2),'--','color',[0.2,0.2,0.2])
end
axis off equal
axis([-0.1,1,0,1])
set(gca,'position',[0 0 1 1],'units','normalized')
print('../../../../master/thesis/fig/ch02/delaunayVoronoiDualityDelVor','-depsc')


