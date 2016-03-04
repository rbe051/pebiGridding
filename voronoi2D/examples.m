% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%% Single fault intersected by several wells 
% close all, clear
% 
% wellLine = {[0.6,0.2;0.65,0.6],...        
%             [0.3,0.3;0.7,0.8],...
%             [0.6,0.2;0.85,0.4],...
%             [0.15,0.7;0.4,0.7]};
% %wellLine ={[0.3,0.3;0.7,0.8]}; 
% 
% x = linspace(0.2, 0.8, 10);
% y = 0.8- 0.5*x - 0.05* sin(6*pi*x);
% fracture = {[x' , y']};
% 
% pri = [2,3,1,4,5];
%    
% Gp = compositePebiGrid(1/24, [1, 1], ...
%                        'wellLines', wellLine, 'wellGridFactor', 0.5^2, ...
%                        'faultLines',fracture, 'faultGridFactor', 1/sqrt(2),...
%                        'circleFactor', 0.6,'mlqtMaxLevel', 2, ...
%                        'mlqtLevelSteps',[0.06,0.025]', 'priOrder', pri);
% 
% Gdist = pebiGrid(1/24, [1, 1], 'wellLines', wellLine, ...
%                 'wellGridFactor', 0.5^2, 'wellRefinement', true, ...
%                 'epsilon',1/12, ...
%                 'faultlines', fracture, 'circleFactor', .6,...
%                 'faultGridFactor', 1/sqrt(2),'priOrder', pri);
%                             
%% Complex wells intersecting
close all; clear
fracture = {};
x = linspace(0.2,0.8);
wellLine = {[0.5,0.2; 0.5,0.3;0.47,0.4;0.4,0.5; 0.33,0.6;0.26,0.7], ...
            [0.5,0.3;0.53,0.4;0.58,0.5],...            
            [0.5,0.45;0.5,0.55;0.45,0.65;0.4,0.75;0.38,0.85],...
            [0.5,0.55;0.55,0.65;0.6,0.75;0.62,0.85]};
                        

Gp = compositePebiGrid(1/14, [1, 1], 'wellLines', wellLine, ...
                      'wellGridFactor', 0.5^2, ...
                      'mlqtMaxLevel', 2, 'mlqtLevelSteps',[0.09,0.04]');
Gdist = pebiGrid(1/14, [1, 1], 'wellLines', wellLine, ...
                'wellGridFactor', 0.5^2, 'wellRefinement',true, 'epsilon',1/7);
                  
                  
              
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
