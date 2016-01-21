clc, clear, close all

% load data
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

% set up data
faultLines = {};
for i = 1:size(x1,1)
    faultLines = {faultLines{:}, [x0(i,:); x1(i,:)]};
end
faultLines = {topCirc, botCirc, recRight,recLeft, faultLines{:}};
wellLines = {circCenter};
% Set parameters
gridSize = 0.5;
xmax = 16;
ymax = 10;
faultGridFactor = 1/sqrt(2); % ratio between grid cells(gridSize) and fault 
                             % grid cells
wellGridFactor = 0.5^4;      % Ratio between grid cells(gridSize) and well 
                             % cells
epsilon = 1;                 % The width of the refinement. Epsilon small 
                             % gives small width

% Generate grid
G = compositeGridPEBIdistmesh(gridSize, [xmax, ymax], 'faultlines', faultLines, ...
                              'faultGridFactor', faultGridFactor, ...
                              'wellLines', wellLines, 'wellRefinement',true,...
                              'epsilon', epsilon,...
                              'wellGridFactor', wellGridFactor);
                          
                          
% plot grid
orange = [1,138/255,0.1];      
figure()
hold on
plotGrid(G, 'facecolor', 'none')
axis equal tight off
hold on
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