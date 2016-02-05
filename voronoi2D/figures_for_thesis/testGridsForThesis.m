addpath('../')

%%
% % curved well pluss two vertical wells
% close all; clear;
% 
% x = linspace(0.2,0.8);
% wellLine = {[0.7105,0.1842], ...
%             [0.2,0.1842], ...
%             [x', 0.5*sin(pi*x)'+0.2]};%, ...
%             %[0.3,0.3;0.7,0.8]}
% fracture = {};
% 
% %With refinement
% % Gp = compositeGridPEBI(1/19, [1, 1], 'wellLines', wellLine, ...
% %                       'wellGridFactor',0.25, 'mlqtMaxLevel', 2, ...
% %                       'mlqtLevelSteps',[0.07,0.03]');
% % Gdist = compositeGridPEBIdistmesh(1/19, [1, 1], 'wellLines', wellLine, ...
% %                              'wellGridFactor', 0.25, 'wellRefDist',1/9);
% %  
% %Without refinement
% Gp = compositeGridPEBI(1/19, [1, 1], 'wellLines', wellLine, ...
%                        'wellGridFactor',0.5, 'mlqtMaxLevel', 0, ...
%                        'mlqtLevelSteps',[0.07,0.03]');
% Gdist = compositeGridPEBIdistmesh(1/19, [1, 1], 'wellLines', wellLine, ...
%                               'wellGridFactor', 0.5, 'wellRefDist',1/500);
%                            
% Gp.cells
% Gdist.cells
% 


%% Single fault intersected by several wells 
close all, clear

wellLine = {[0.6,0.2;0.65,0.6],...        
            [0.3,0.3;0.7,0.8],...
            [0.6,0.2;0.85,0.4],...
            [0.15,0.7;0.4,0.7]};
%wellLine ={[0.3,0.3;0.7,0.8]}; 

x = linspace(0.2, 0.8, 10);
y = 1- 0.8*x - 0.1* sin(6*pi*x);
fracture = {[x' , y']};

pri = [2,3,1,4,5];

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
Gp = compositeGridPEBI(1/24, [1, 1], ...
                       'wellLines', wellLine, 'wellGridFactor', 0.5^2, ...
                       'faultLines',fracture, 'faultGridFactor', 1/sqrt(2),...
                        'circleFactor', 0.6,'mlqtMaxLevel', 2, ...
                        'mlqtLevelSteps',[0.06,0.025]', 'priOrder', pri);

Gdist = compositeGridPEBIdistmesh(1/24, [1, 1], 'wellLines', wellLine, ...
                                'wellGridFactor', 0.5^2, 'wellRefinement', true, ...
                                'epsilon',1/12, ...
                                'faultlines', fracture, 'circleFactor', .6,...
                                'faultGridFactor', 1/sqrt(2),'priOrder', pri);

Gp.cells
Gdist.cells



%% Complex wells intersecting
% close all; clear
% fracture = {};
% x = linspace(0.2,0.8);
% wellLine = {[0.5,0.2; 0.5,0.3;0.47,0.4;0.4,0.5; 0.33,0.6;0.26,0.7], ...
%             [0.5,0.3;0.53,0.4;0.58,0.5],...            
%             [0.5,0.45;0.5,0.55;0.45,0.65;0.4,0.75;0.38,0.85],...
%             [0.5,0.55;0.55,0.65;0.6,0.75;0.62,0.85]};
%                         
% 
% Gp = compositeGridPEBI(1/14, [1, 1], 'wellLines', wellLine, ...
%                       'wellGridFactor', 0.5^2, ...
%                       'mlqtMaxLevel', 2, 'mlqtLevelSteps',[0.09,0.04]');
% Gdist = compositeGridPEBIdistmesh(1/14, [1, 1], 'wellLines', wellLine, ...
%                                  'wellGridFactor', 0.5^2, 'wellRefDist',1/7);
%                   
%                   
% Gp.cells
% Gdist.cells
                            
%% Plotting                       
orange = [1,138/255,0.1];      
figure()
hold on
plotGrid(Gp, 'faceColor', 'none')
axis equal tight off
hold on
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
%plotFault(Gp)
plotWells(Gp)
% % close up
% box = [0.55,0.75,0.25,0.45];
% rectangle('position', [box(1),box(3),box(2)-box(1),box(4)-box(3)], 'linewidth', 2)

figure()
hold on
plotGrid(Gdist, 'faceColor', 'none')
axis equal tight off
hold on
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
%plotFault(Gdist)
plotWells(Gdist)
%close up
%rectangle('position', [box(1),box(3),box(2)-box(1),box(4)-box(3)], 'linewidth', 2)




