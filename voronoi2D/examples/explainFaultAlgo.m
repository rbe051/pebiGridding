% This script goes trough the most important details about how 
% createFaultPoints create the fault points.
% It does not cover the details, nor does it cover intersecting faults, but
% it will hopefully enlighten persons that whish to implement their own
% algorithms conforming to faults. See explainWellAlgo for creating well
% points.


% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.

%clear all; close all


%% Create a fault
fault    = [0.2,0.8; 0.5,0.5; 0.8,0.5];  % A fault with a single corner
dFault   = 0.1;                          % distance between circles
circFact = 0.6;                          % a constant that dicides the 
                                         % radii of the circles
%% Plot Fault
figure()
hold on
plot(fault(:,1), fault(:,2));
axis equal tight
axis([0.1,0.9,0.4,0.9])

%% Interpolate fault
% We interpolate the fault with a set of equidistant points. We set the
% circlecenters to these points
linesDist = sqrt(sum(diff(fault,[],1).^2,2));  % Total length of each line segment
linesDist = [0; linesDist];                    % Add the starting point
cumDist   = cumsum(linesDist);                 % find cummulative distant to each point
dt        = cumDist(end)/...                   % Make sure end point is included
              ceil(cumDist(end)/dFault);       
newPtsT   = 0:dt:cumDist(end);
CC        = interp1(cumDist, fault, newPtsT);  % Circle centers

%% Plot interpolation:
yellow =  [0.9290    0.6940    0.1250];
plot(CC(:,1),CC(:,2),'.','color',yellow,'markersize',15);


%% Calculate circle radiuses
d  = sqrt(sum((CC(2:end,:)-CC(1:end-1,:)).^2, 2)); % Find the distance between all circles
CR = circFact*[d(1); (d(1:end-1) + d(2:end))/2; d(end)];
                                                   % Set radius to the
                                                   % average of the line
                                                   % segments on both sides
                                                   % of the circle
            
            
%% Calculate midpoints and circle offsets
a = (d.^2 - CR(2:end).^2 + CR(1:end-1).^2)./(2*d);
h = sqrt(CR(1:end-1).^2 - a.^2);
    
% Calculate normal vectors
n1 = (CC(2:end,:)-CC(1:end-1,:))./repmat(d,1,2); %Unit vector
n2 = [-n1(:, 2), n1(:,1)];                       %Unit normal


%% Calculate intersections
lPts = CC(1:end-1,:) + bsxfun(@times, a, n1)  ...
                     + bsxfun(@times, h,  n2);
rPts = CC(1:end-1,:) + bsxfun(@times, a, n1)  ...
                     - bsxfun(@times, h,  n2);
pts = [lPts; rPts];

%% Plotting
yellow =  [0.9290    0.6940    0.1250];
blue   = [0    0.4470    0.7410];
theta= linspace(0,2*pi)';
for i = 1:size(CC,1)
  X = repmat(CC(i,:),100,1) + repmat(CR(i),100,2).*[cos(theta), sin(theta)];
  plot(X(:,1),X(:,2),'k')
end
plot(pts(:,1),pts(:,2),'.','color',blue, 'markersize',15)
%% Cartesian background-grid
[X,Y]   = meshgrid(0:dFault:1);
resPts = [X(:),Y(:)];
CRrep   = repmat(CR',size(resPts,1),1);
removed = any(pdist2(resPts,CC)<CRrep,2); % Remove any resPts innside one of the 
resPts = resPts(~removed,:);              % circles calculated above

%% Plot all seeds
plot(resPts(:,1), resPts(:,2),'.k','markersize',15);


%% Create full grid
Gt = triangleGrid([pts;resPts]);
G  = pebi(Gt);
figure()
plotGrid(G,'facecolor','none')
axis equal tight
axis([0.1,0.9,0.4,0.9])









