clc; clear all; close all

% set boundary
boundary = [-1,-1,-1;  ...
             1,-1,-1;  ...
             1, 1,-1;  ...
            -1, 1,-1;  ...
            -1,-1, 1;  ...
             1,-1, 1;  ...
             1, 1, 1;  ...
            -1, 1, 1]; ...
%             
% boundary = [-1, 0, -1/sqrt(2);...
%              1, 0, -1/sqrt(2);...
%              0,-1,  1/sqrt(2);...
%              0, 1,  1/sqrt(2)]; ...
%              
face = [-4,-4;-4,1;1,-4;1,1]*1.1;

face = delaunayTriangulation(face);
faceDt.Points = [face.Points,zeros(size(face.Points,1),1)];
%faceDt.Points = [face.Points, [-0.2;0.2;0.2;-0.2]];
faceDt.ConnectivityList = face.ConnectivityList;

% initilize random points
n = 200;

pts = rand(n,3);
pts(:,1) = pts(:,1)*(max(boundary(:,1))-min(boundary(:,1))) + min(boundary(:,1));
pts(:,2) = pts(:,2)*(max(boundary(:,2))-min(boundary(:,2))) + min(boundary(:,2));
pts(:,3) = pts(:,3)*(max(boundary(:,3))-min(boundary(:,3))) + min(boundary(:,3));

% pts = pts + 0.0*rand(size(pts));
% [G,optPts,f,g] = optiVoronoi3D(pts,boundary);
% plotGrid(G)
% figure()
  boundary = delaunayTriangulation(boundary);
  
  G = restrictedVoronoiDiagram(pts,boundary);
  figure()
  plotGrid(G);
  
  [Gc,optPtsc,fc,gc] = optiVoronoi3DClip(pts,boundary,'fault', faceDt,...
                                         'tol', 1e-6);
  figure()
  plotGrid(Gc)
  
%save('unitSquare200')
    
