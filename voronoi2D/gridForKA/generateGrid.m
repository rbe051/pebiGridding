clc, clear, close all

%% load data
load('data.mat')
% convert to cartesian
x1 = x0 + polarToCart(theta, r);

% Plot data
figure
hold on
for i = 1:size(x0,1)
    plot([x0(i,1), x1(i,1)], [x0(i,2), x1(i,2)], 'k')
end
plot(recRight(:,1), recRight(:,2));
plot(recLeft(:,1), recLeft(:,2));
plot(topCirc(:,1), topCirc(:,2));
plot(botCirc(:,1), botCirc(:,2));
axis equal
figure()

%% set up data
faultLines = {};
for i = 1:size(x1,1)
    faultLines = {faultLines{:}, [x0(i,:); x1(i,:)]};
end
% The well is set as a fault to make edges in the grid follow its boundary
faultLines = {topCirc, botCirc, recRight,recLeft, faultLines{:}};
% This is only needed if you want refinement around the well
wellLines = {circCenter};

% Set parameters
gridSize = 0.8;
xmax = 16;
ymax = 10;
faultGridFactor = 1/sqrt(2); % ratio between grid cells(gridSize) and fault 
                             % grid cells
wellGridFactor = 0.5^4;      % Ratio between grid cells(gridSize) and well 
                             % cells
epsilon = 1;                 % The width of the refinement. Epsilon small 
                             % gives small width

%% Generate grid
G = pebiGrid(gridSize, [xmax, ymax], 'faultlines', faultLines, ...
            'faultGridFactor', faultGridFactor, 'wellLines', wellLines, ...
            'wellRefinement',true, 'epsilon', epsilon, ...
            'wellGridFactor', wellGridFactor);
                     
% Fix labeling of cells and edges. The well and perforation cells is labled
% as G.cells.type = true. If you wish to treat them differently you can use
% the logical indexes isWell and isPerf
G = computeGeometry(G);
centr = G.cells.centroids;
% Lable cells inside well
isWell = (circCenter(1)- centr(:,1)).^2 + (circCenter(2) - centr(:,2)).^2 < circR^2;
G.cells.tag(isWell) = true;
% labele cells inside rectangle
inLeftRec  = inpolygon(centr(:,1), centr(:,2), recLeft(:,1), recLeft(:,2))...
            & (circCenter(1)- centr(:,1)).^2 + (circCenter(2) - centr(:,2)).^2 > circR^2;
inRightRec = inpolygon(centr(:,1), centr(:,2), recRight(:,1), recRight(:,2))...
            & (circCenter(1)- centr(:,1)).^2 + (circCenter(2) - centr(:,2)).^2 > circR^2;
isPerf = inLeftRec | inRightRec;
G.cells.tag(isPerf) = true;

% Remove unwanted fault tags.
indx = find(isPerf | isWell);
remTag = false(G.faces.num);
for i = 1:size(indx,1)
    facePos = G.cells.facePos(indx(i)):G.cells.facePos(indx(i)+1)-1;
    remTag(G.cells.faces(facePos)) = true;
end
G.faces.tag(remTag) = false;

%% plot grid
orange = [1,138/255,0.1];      
figure()
hold on
plotGrid(G, 'facecolor', 'none')
axis equal tight off
hold on
%plot(G.cells.centroids(inRec,1), G.cells.centroids(inRec,2), '.')
for i = 1:numel(wellLines)
  line = wellLines{i};
  if size(line,1) == 1
      plot(line(1,1), line(1,2),'.r', 'markersize', 8);
  end
  plot(line(:, 1), line(:, 2),'r');
end
for i = 1:numel(faultLines)
  line = faultLines{i};
  plot(line(:, 1), line(:, 2),'--', 'color',orange);
end