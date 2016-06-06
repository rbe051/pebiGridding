clear; close all
%addpath ../
%addpath ../../../vem/mat/VEM2D
%addpath ../../../vem/mat/
%addpath ../../voronoi3D/
%% Set boundary and parameters
k = 100;
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

G = clippedPebi2D(pts, x0);
G = computeGeometry(G);

%% Plotting
col = get(gca,'ColorOrder');
box = [-0.1,0;-0.1,1;1,1;1,0;-0.1,0];

%% Only voronoi
figure(); hold on
plotGrid(G,'facecolor','none')
plot(pts(:,1), pts(:,2),'.','color',col(1,:),'markersize',20);
circCenters = G.nodes.coords;
plot(circCenters(:,1), circCenters(:,2),'.','color',col(4,:),'markersize',20);
plot(box(:,1), box(:,2),'k--');
axis off equal

axis([-0.1,1,0,1])
%set(gca,'position',[0 0 1 1],'units','normalized')
print('../../../../master/thesis/fig/ch02/delaunayVoronoiDualityVor','-depsc')

%% Only delaunay
figure(); hold on
plotGrid(Gt,'facecolor','none')
plot(pts(:,1), pts(:,2),'.','color',col(1,:),'markersize',15);
plot(circCenters(:,1), circCenters(:,2),'.','color',col(4,:),'markersize',20);
% Plot circumcircle
t = 5;
dt = delaunayTriangulation(pts);
[CC, r] = circumcenter(dt,t);
theta = linspace(0,2*pi)';
for i = 1:size(CC)
X = repmat(CC(i,:),100,1) + repmat(r(i),100,2).*[cos(theta), sin(theta)];
plot(X(:,1), X(:,2),'--','color',[0.2,0.2,0.2])
end
plot(box(:,1), box(:,2),'k--');
axis off equal
axis([-0.1,1,0,1])
%set(gca,'position',[0 0 1 1],'units','normalized')
print('../../../../master/thesis/fig/ch02/delaunayVoronoiDualityDel','-depsc')


%% Both delaunay and voronoi
figure(); hold on
plotGrid(Gt,'facecolor','none')
plotGrid(G,'facecolor','none')
plot(pts(:,1), pts(:,2),'.','color',col(1,:),'markersize',20);
plot(circCenters(:,1), circCenters(:,2),'.','color',col(4,:),'markersize',20);
% Plot circumcircle
t = 5;
dt = delaunayTriangulation(pts);
[CC, r] = circumcenter(dt,t);
theta = linspace(0,2*pi)';
for i = 1:size(CC)
X = repmat(CC(i,:),100,1) + repmat(r(i),100,2).*[cos(theta), sin(theta)];
plot(X(:,1), X(:,2),'--','color',[0.2,0.2,0.2])
end
plot(box(:,1), box(:,2),'k--');
axis off equal
axis([-0.1,1,0,1])
%set(gca,'position',[-0 0 1 1],'units','normalized')
print('../../../../master/thesis/fig/ch02/delaunayVoronoiDualityDelVor','-depsc')


