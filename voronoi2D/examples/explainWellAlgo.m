% This script goes trough the most important details about how 
% createWellPoints works.
% It does not cover the details, nor does it cover intersecting wells, but
% it will hopefully enlighten persons that whish to implement their own
% algorithms conforming to wells. 
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.


%% Set path to voronoi2d
%add path ../
%clear; close all

%% Create a Well
% We create a circular well to demonstrate the flexibillity
% of the gridding method. 
dWell = 0.1;                         % Grid cell size
n     = 20;                          % Number of line-segments that approximate the well
dT    = 2*pi/20;
theta = 2*pi*(0:1/n:1 -1/n)';       
well  = repmat([0.5,0.5],n,1) ...
          + 0.3*[cos(theta), sin(theta)];

%% Plot well
figure()
hold on
plot([well(:,1);well(1,1)], [well(:,2);well(1,2)]);
axis equal tight
axis([0.1,0.9,0.1,0.9])
%% Interpolate well
% We interpolate the well with a set of equidistant points
linesDist = sqrt(sum(diff(well,[],1).^2,2)); % Total length of each line segment
linesDist = [0; linesDist];                  % Add starting points
cumDist   = cumsum(linesDist);               % Find cummulative distant to each segment
dt        = cumDist(end)/...                 % Make sure the end point is included
              ceil(cumDist(end)/dWell);
newPtsT   = 0:dt:cumDist(end);               % interpolation points
wPts   = interp1(cumDist, well, newPtsT);    % well seeds

%% Plot interpolation:
yellow = [0.9290    0.6940    0.1250];
blue   = [0    0.4470    0.7410];
plot(wPts(:,1),wPts(:,2),'.','color',blue,'markersize',15);

%% Cartesian background-grid
[X,Y]   = meshgrid(0:dWell:1);
resPts = [X(:),Y(:)];
removed = any(pdist2(resPts,wPts)<dWell,2);  % Remove any resPts that are 
resPts = resPts(~removed,:);                 % to close to a well seed

plot(resPts(:,1), resPts(:,2),'.k','markersize',15);

%% Create full grid
Gt = triangleGrid([wPts;resPts]);
G  = pebi(Gt);
figure()
plotGrid(G,'facecolor','none')
axis equal tight
axis([0.1,0.9,0.1,0.9])






