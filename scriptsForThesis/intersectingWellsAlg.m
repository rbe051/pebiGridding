clc; clear;close all
addpath ../voronoi2D
blue   = [0       0.4470, 0.7410];
red    = [0.8500, 0.3250, 0.0980];
yellow = [0.9290,  0.6940, 0.1250];

wells = {[0.5,0.2;0.5,0.8],[0.2,0.3;0.8,0.8]};
ds = 0.1;

[wells] = splitAtInt(wells, {});

pts = createWellGridPoints(wells,ds);
[X,Y] = meshgrid(0:ds:1,0:ds:1);
resPts = [X(:),Y(:)];

resPts = removeConflictPoints2(resPts, pts, 0.1*ones(size(pts,1),1));

Gt = triangleGrid([pts;resPts]);
G  = pebi(Gt);

w1 = logical([1;1;1;1;0;0;0;0;0;0;0;0;1;1;1;1]);
w2 = logical([0;0;0;0;1;1;1;1;0;1;1;1;0;0;0;0]);
ws = ~w1 & ~w2;

figure()
hold on
for i = 1:numel(wells)
  w = wells{i};
  plot(w(:,1),w(:,2),'k')
end
  
plot(pts(w1,1),pts(w1,2),'.','color',blue,'markersize',20);
plot(pts(w2,1),pts(w2,2),'.','color',red,'markersize',20);
plot(pts(ws,1),pts(ws,2),'.','color',yellow,'markersize',20);

axis([0.3,0.7,0.35,0.75])
axis off

figure(); hold on
plotGrid(G,'facecolor','none')
plot(pts(w1,1),pts(w1,2),'.','color',blue,'markersize',20);
plot(pts(w2,1),pts(w2,2),'.','color',red,'markersize',20);
plot(pts(ws,1),pts(ws,2),'.','color',yellow,'markersize',20);

axis([0.3,0.7,0.35,0.75])
axis off
