%% Single fault intersected by several wells 
close all

wellLine = {[0.5,0.2; 0.5,0.35;0.47,0.6;0.4,0.75; 0.33,0.85;0.26,0.9], ...
            [0.8,0.15; 0.3,0.5]};
        
fracture = {[0.2,0.8;0.8,0.2]};

%% Plott wells and Fault
figure()
orange = [1,138/255,0.1];    
hold on
axis([0,1,0,1])
axis off
for i = 1:numel(wellLine)
  line = wellLine{i};
  if size(line,1) == 1
      plot(line(1,1), line(1,2),'.r', 'markersize', 8);
  end
  plot(line(:, 1), line(:, 2),'r');
end
for i = 1:numel(fracture)
  line = fracture{i};
  plot(line(:, 1), line(:, 2),'color', orange);
end



% Without refinement
% Gp = compositeGridPEBI(1/24, [1, 1], ...
%                        'wellLines', wellLine, 'wellGridFactor', 24/26/2, ...
%                        'faultLines',fracture, 'faultGridFactor', 1/sqrt(2),...
%                         'circleFactor', 0.6);
% Gdist = compositeGridPEBIdistmesh(1/24, [1, 1], 'wellLines', wellLine, ...
%                                 'wellGridFactor', 0.5*24/26, 'wellRefDist',1/500, ...
%                                 'faultlines', fracture, 'circleFactor', .6,...
%                                 'faultGridFactor', 1/sqrt(2));

% Whith refinement      
%Gp = compositeGridPEBI(1/16, [1, 1], ...
%                       'wellLines', wellLine, 'wellGridFactor', 0.5^2, ...
%                       'faultLines',fracture, 'faultGridFactor', 1/sqrt(2),...
%                        'circleFactor', 0.6,'mlqtMaxLevel', 2, ...
%                        'mlqtLevelSteps',[0.075,0.035]',...
%                        'fullFaultEdge', 1);


Gdist = compositeGridPEBIdistmesh(1/16, [1, 1], 'wellLines', wellLine, ...
                                'wellGridFactor', 0.5^2, 'wellRefDist',1/12, ...
                                'faultlines', fracture, 'circleFactor', .6,...
                                'faultGridFactor', 1/sqrt(2));

Gp.cells
Gdist.cells



%% Plotting                       
orange = [1,138/255,0.1];      
figure()
hold on
plotGrid(Gp, 'faceColor', 'none')
axis equal tight off
hold on
%plotFault(Gp)
plotWells(Gp)
for i = 1:numel(wellLine)
  line = wellLine{i};
  if size(line,1) == 1
      plot(line(1,1), line(1,2),'.r', 'markersize', 8);
  end
  plot(line(:, 1), line(:, 2),'r');
end
for i = 1:numel(fracture)
  line = fracture{i};
  plot(line(:, 1), line(:, 2),'color',orange);
end


figure()
hold on
plotGrid(Gdist, 'faceColor', 'none')
axis equal tight off
hold on
%plotFault(Gdist)
plotWells(Gdist)
for i = 1:numel(wellLine)
  line = wellLine{i};
  if size(line,1) == 1
      plot(line(1,1), line(1,2),'.r', 'markersize', 8);
  end
  plot(line(:, 1), line(:, 2),'r');
end
for i = 1:numel(fracture)
  line = fracture{i};
  plot(line(:, 1), line(:, 2),'color', orange);
end
