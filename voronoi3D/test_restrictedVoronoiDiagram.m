clc; clear; close all


% %% Cube
% dtPts = [-1,-1,-1;  ...
%           1,-1,-1;  ...
%           1, 1,-1;  ...
%          -1, 1,-1;  ...
%          -1,-1, 1;  ...
%           1,-1, 1;  ...
%           1, 1, 1;  ...
%          -1, 1, 1]; ...
% %dtPts = [dtPts; 2*rand(4,3)-1];
% 
% 
% pts = 2*rand(20,3)-1;
% % pts=    [-0.5443    0.4773    0.2519;...
% %          -0.0038    0.1720    0.3219;...
% %           0.8017   -0.5065    0.4595;...
% %           0.1493    0.3328    0.7815;...
% %           0.6904   -0.8330    0.9646];...
% %% ball
n = 10;
[dtx, dty, dtz] = sphere(30);
dtPts = [dtx(:), dty(:), dtz(:)];
[ptsx,ptsy, ptsz] = sphere(n);
pts = 0.5*[ptsx(:), ptsy(:), ptsz(:)];
pts = pts + 0.2*rand((n+1)^2,3);

dt = delaunayTriangulation(dtPts);
%% sinus
% n = 10;
% r1 = 1;
% r2 = 2;
% xb = linspace(-1,1,10);
% 
% [xb,yb,zb] = ndgrid(xb,xb,[1,-1]);
% dtpts = [xb(:),yb(:),zb(:)];
% 
% x  = linspace(-1,1,n);
% [x,y,z] = ndgrid(x,x,0);
% pts = 0.8*[x(:),y(:),z(:)];
% pts = unique(pts,'rows');
% pts = pts + 0.4*rand(size(pts));
% pts(:,3) = pts(:,3) + 0.2*sin(pi*pts(:,1));
% 
% dt = delaunayTriangulation(dtpts);
% 
% for i = 1:size(dt.ConnectivityList,1)
%     hull = convhulln(dt.Points(dt.ConnectivityList(i,:),:));
%     patch('Vertices', dt.Points(dt.ConnectivityList(i,:),:), 'faces', hull,'facecolor','y')
% end
% figure()
% 
% dt = struct('Points', dt.Points,...
%             'ConnectivityList',dt.ConnectivityList,...
%             'freeBoundary', dt.freeBoundary,...
%             'edges', dt.edges);
% dt.Points(:,3) = dt.Points(:,3) + 0.5*sin(pi*dt.Points(:,1));

%%
% for i = 1:size(dt.ConnectivityList,1)
%     hull = convhulln(dt.Points(dt.ConnectivityList(i,:),:));
%     patch('Vertices', dt.Points(dt.ConnectivityList(i,:),:), 'faces', hull,'facecolor','y')
% end
t = tic
G = restrictedVoronoiDiagram(pts, dt);
toc(t)

% for i = 1:G.cells.num
%     plotGrid(G, i, 'edgecolor','r','facealpha',0.2);
% end
