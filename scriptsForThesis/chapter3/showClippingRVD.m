clear; close all




%% Creat domain
dtPts = [-1,-1,-1;  ...
          1,-1,-1;  ...
          1, 1,-1;  ...
         -1, 1,-1;  ...
         -1,-1, 1;  ...
          1,-1, 1;  ...
          1, 1, 1;  ...
         -1, 1, 1];

bound = delaunayTriangulation(dtPts);
%% Create sites
n = 5;
dx = 0.1;
x = linspace(-1+dx,1-dx,n);
[X,Y,Z] = ndgrid(x,x,x);

X(:,1:2:end) = X(:,1:2:end) + 1/n/sqrt(2);
Y(:,2:2:end) = Y(:,2:2:end) + 1/n/sqrt(2);
pts = [X(:), Y(:), Z(:)];
pts = pts;

G = restrictedVoronoiDiagram(pts, bound);

plotGrid(G)

x = linspace(-1+dx,1-dx,n);
[X,Y] = meshgrid(x,x);

ptri = [X(:),Y(:)];
ptri = ptri + 0.1*rand(size(ptri));
tri = delaunayTriangulation(ptri);

patch('Vertices', [tri.Points,-1.1*ones(size(tri.Points,1),1)], 'faces', tri.ConnectivityList,'facealpha',0.2, ...
      'facecolor', 'none')