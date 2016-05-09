clear; close all
addpath ../../../vem/mat/VEM2D
addpath ../../../vem/mat/

%% Set boundary and parameters
k = 100;
n = [-1, 0; ...
      0, 1; ...
      1, 0; ...
      0,-1];
x0 = [0, 0; ...
      0, 1; ...
      1, 1; ...
      1, 0];
    

%% Create none CVD
pts = rand(k,2);


G = clippedPebi2D(pts, n, x0);
G = computeGeometry(G);

%% Create CVD
% row1 = [linspace(0.165,0.835,3)',repmat(0.165,3,1)];
% row2 = [linspace(.12,.88,4)',repmat(0.5,4,1)];
% row3 = [linspace(0.165,0.835,3)',repmat(0.835,3,1)];
% 
% pCVT = [row1;row2;row3];
% Gcvt = clippedPebi2D(pCVT,n,x0);
% Gcvt = computeGeometry(Gcvt);

[Gcvt,pCVT] = createCVD(k,n,x0);
Gcvt = computeGeometry(Gcvt);
%% Plotting
figure(); hold on
plotGrid(G,'facecolor','none')
plot(pts(:,1), pts(:,2),'.k','markersize',8);
plot(G.cells.centroids(:,1), G.cells.centroids(:,2),'ok','markersize',8)
axis off equal tight
%print('../../../../master/thesis/fig/ch02/nonCVD','-deps')

figure(); hold on
plotGrid(Gcvt,'facecolor','none')
plot(pCVT(:,1), pCVT(:,2),'.k','markersize',8);
plot(Gcvt.cells.centroids(:,1), Gcvt.cells.centroids(:,2),'ok','markersize',8)
axis off equal tight

%print('../../../../master/thesis/fig/ch02/CVD','-deps')

