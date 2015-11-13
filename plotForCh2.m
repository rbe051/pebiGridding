clear all; close all
pts = [0.65,-0.2;-0.5,0.7;0,0; 0.3, 0.7; -0.2,1.25; 0.85,1.8; 1.1,0.5; 0.5,1.2];
%pts = 10*randn(6,2);
x = 1:3;
y = 1:3;
[X,Y] = meshgrid(x,y);
int = 1:size(X,1);
X(int,int) = X(int,int) + 1*randn(size(X(int,int)));
Y(int,int) = Y(int,int) +1*randn(size(X(int,int)));
%X(2,1) = 0;
%pts = [X(:),Y(:)];
G =delaunayTriangulation(pts);
Gt = triangleGrid(G.Points, G.ConnectivityList)

% Add points at infinite before creating the pebi grid
pts = [pts;-100,-100;-100,100;100,100;100,-100];
Gp =delaunayTriangulation(pts);
Gp = triangleGrid(Gp.Points, Gp.ConnectivityList);
Gp = pebi(Gp);
[CC,r] = circumcenter(G);

circles = {};
n = 100;
theta = linspace(0,2*pi,n)';
for i = 1:length(CC)
   circles{i} = repmat(CC(i,:),n,1) + r(i) * [cos(theta),sin(theta)];
end
axis_lim = [-0.5,1.2,-0.5,2];
plotGrid(Gt, 'facecolor', 'none')
hold on
plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
axis(axis_lim)
axis equal off
figure()
hold on
plotGrid(Gp,'facecolor','none')
plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
axis(axis_lim)
axis equal off
figure()
plotGrid(Gt, 'facecolor', 'none')
axis(axis_lim)
hold on
plotGrid(Gp,'facecolor','none')
plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
axis(axis_lim)
for i = 1:numel(circles)
   y = circles{i}
   %plot(y(:,1), y(:,2), 'color', [0.7,0.7,0.7])
end

%plot(CC(:,1), CC(:,2),'o')
%plot(Gp.nodes.coords(:,1), Gp.nodes.coords(:,2),'ro')
axis equal off