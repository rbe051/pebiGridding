%% Example of duality of delaunay and Voronoi.
% clear all; close all
% pts = [0.65,-0.2;-0.5,0.7;0,0; 0.3, 0.7; -0.2,1.25; 0.85,1.8; 1.1,0.5; 0.5,1.2];
% %pts = 10*randn(6,2);
% x = 1:3;
% y = 1:3;
% [X,Y] = meshgrid(x,y);
% int = 1:size(X,1);
% X(int,int) = X(int,int) + 1*randn(size(X(int,int)));
% Y(int,int) = Y(int,int) +1*randn(size(X(int,int)));
% %X(2,1) = 0;
% %pts = [X(:),Y(:)];
% G =delaunayTriangulation(pts);
% Gt = triangleGrid(G.Points, G.ConnectivityList)
% 
% % Add points at infinite before creating the pebi grid
% pts = [pts;-100,-100;-100,100;100,100;100,-100];
% Gp =delaunayTriangulation(pts);
% Gp = triangleGrid(Gp.Points, Gp.ConnectivityList);
% Gp = pebi(Gp);
% [CC,r] = circumcenter(G);
% 
% circles = {};
% n = 100;
% theta = linspace(0,2*pi,n)';
% for i = 1:length(CC)
%    circles{i} = repmat(CC(i,:),n,1) + r(i) * [cos(theta),sin(theta)];
% end
% axis_lim = [-0.5,1.2,-0.5,2];
% plotGrid(Gt, 'facecolor', 'none')
% hold on
% plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
% axis(axis_lim)
% axis equal off
% figure()
% hold on
% plotGrid(Gp,'facecolor','none')
% plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
% axis(axis_lim)
% axis equal off
% figure()
% plotGrid(Gt, 'facecolor', 'none')
% axis(axis_lim)
% hold on
% plotGrid(Gp,'facecolor','none')
% plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
% axis(axis_lim)
% for i = 1:numel(circles)
%    y = circles{i}
%    %plot(y(:,1), y(:,2), 'color', [0.7,0.7,0.7])
% end
% 
% %plot(CC(:,1), CC(:,2),'o')
% %plot(Gp.nodes.coords(:,1), Gp.nodes.coords(:,2),'ro')
% axis equal off


%% Showing not sufficient condition
% clear; close all
% faultLine = {[1,0.3;1,1.7]};
% orange = [1,138/255,0.1];      
% 
% Gp = compositeGridPEBI(0.333,[2,2], 'faultLines', faultLine, ...
%                        'faultGridFactor', 2/0.667,'circleFactor',0.55,...
%                        'fullFaultEdge', 0)
% plotGrid(Gp,'facecolor','none')
% line = faultLine{1};
% plot(line(:,1), line(:,2), 'color', orange); 
% axis off equal
% Gp = compositeGridPEBI(0.333,[2,2], 'faultLines', faultLine, ...
%                        'faultGridFactor', 2/0.667,'circleFactor',0.55,...
%                        'fullFaultEdge', 1)
% plotGrid(Gp, 'facecolor','none')
% plot(line(:,1), line(:,2), 'color', orange); 
% axis off equal


%% Multilevel quad tree refinement.
close all; clear
dx = 0.1;
x = 0:dx:1;
[X, Y] = meshgrid(x,x);
pts = [X(:), Y(:)];


centerRef = [0.5,0.5];
levels = [0.4,0.2,0.1]';
res = {};
varArg = {'level', 1, 'maxLev', 3, 'distTol', levels};
for i = 1:size(pts,1)
    res = [res; mlqt(pts(i,:), centerRef, dx, varArg{:})];
end
pts = vec2mat([res{:,1}],2);

Tri = delaunayTriangulation(pts);
Gt = triangleGrid(pts, Tri.ConnectivityList);
Gp = pebi(Gt);

Gp = pebi(Gt);
plotGrid(Gp,'facecolor','none')
axis equal off








