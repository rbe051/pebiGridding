clear; close all
% addpath ../../../vem/mat/VEM2D
% addpath ../../../vem/mat/
% addpath ../

%% Set boundary and parameters
k = 200;
n = [-1, 0; ...
      0, 1; ...
      1, 0; ...
      0,-1];
x0 = [0, 0; ...
      0, 1; ...
      1, 1; ...
      1, 0];
    

%% Create CVD
[G,p] = createCVD(k,n,x0);
G = computeGeometry(G);
%% Plotting


figure(); hold on
plotGrid(G,'facecolor','none')

axis off equal tight
