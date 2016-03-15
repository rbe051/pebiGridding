% Multiple fractures. From Fung et.al 15.
% close all
% x = 0.2:0.05:0.8;
% l = {[0.65,0.1;0.65,0.926],...
%      [0.2,0.175; 0.875,0.875], ...
%      [0.2,0.925; 0.9,0.125], ...
%      [0.45,0.15; 0.83, 0.35]};
% 
% Gp = compositeGridPEBI([0.1,-1,-1], [1, 1], 'faultLines', l, 'padding', 1);
% 
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight
% hold on
% plotFault(Gp)
% for i = 1:numel(l)
%   line = l{i};
%   plot(line(:, 1), line(:, 2),'r');
% end
% 


% %close all
% x = 0.2:0.05:0.8;
% l = {[0.65,0.1;0.65,0.926],...
%      [0.2,0.175; 0.875,0.875], ...
%      [0.2,0.925; 0.9,0.125], ...
%      [0.45,0.15; 0.83, 0.35]};
% 
% Gp = compositeGridPEBI([30,30], [1, 1], 'fracLines', l, 'padding', 1);
% 
% 
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight
% hold on
% plotFault(Gp)
% for i = 1:numel(l)
%   line = l{i};
%   %plot(line(:, 1), line(:, 2),'r');
% end


%%
%%Close up on one iregular fracture.
% close all; clear all
% x = [0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8];
% y = [0.25,0.26,0.28,0.34,0.39,0.4,0.39,0.36];
% l = {[x',y']};    
% 
% Gp = compositeGridPEBI([0.1,-1,-1], [1, 0.75], 'faultLines', l, 'padding', 0,...
%                        'fullFaultEdge', 0, 'circleFactor', 0.6);
% 
% figure
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% plotFault(Gp)


% %% Showing sufficient condition for fault edge with circle factor 0.6
% close all
% x = [3, 6, 9, 12, 15];
% y = [7.075, 8, 8.925, 8, 7.075];
% l = {[x',y']};    
% 
% Gp = compositeGridPEBI([9,7], [18, 15.74], 'fracLines', l, 'padding', 1,...
%                        'fullFaultEdge', 0, 'circleFactor', 0.6, ...
%                        'fracGridSize', 4);
% 
% 
% plotGrid(Gp, 'faceColor', 'none')
% 
% axis  equal tight off
% hold on
% plotFault(Gp)
% for i = 1:numel(l)
%   line = l{i};
%   %plot(line(:, 1), line(:, 2),'r');
% end

% %%Complex faults intersecting.
%close all

x5 = linspace(0.35,0.2,10);
y5 = [0.2,0.25,0.3,0.34,0.40,0.5,0.6,0.65,0.7,0.85];
x4 = linspace(0.2,0.7,10);
y4 = [0.45,0.4,0.4,0.38,0.35,0.35,0.4,0.5,0.6,0.73];
x3 = linspace(0.1,0.9,10);
y3 = [0.8,0.85,0.85,0.9,0.9,0.85,0.8,0.70,0.7,0.8];
x2 = linspace(0.8,0.55,10);
y2 = [0.2,0.25,0.3,0.34,0.40,0.5,0.6,0.65,0.7,0.85];
x1 = linspace(0.2,0.7,10);
y1 = [0.6,0.57,0.56,0.52,0.50,0.5,0.53,0.57,0.59,0.63];


l = {[x2',y2'], [x3',y3'],  [x5',y5'],[x1',y1'],[x4',y4']};    

%l = {[0.4,0.2;0.8,0.8],[0.5,0.25;0.6,0.9],[0.4,0.7;0.8,0.65]};    

Gp = compositePebiGrid(1.155/50, [1, 1], 'faultLines', l, 'faultGridFactor', 1/1.155,...
                        'fullFaultEdge', 1, 'circleFactor', 0.6);

plotGrid(Gp, 'faceColor', 'none')
axis equal tight off
hold on
%plotFault(Gp)

for i = 1:numel(l)
  line = l{i};
  plot(line(:, 1), line(:, 2),'r');
end

% % Multiple fractures. From Fung et.al 15. huge frac gird size
% close all
% x = 0.2:0.05:0.8;
% l = {[0.65,0.1;0.65,0.926],...
%      [0.2,0.175; 0.875,0.875], ...
%      [0.2,0.925; 0.9,0.125], ...
%      [0.45,0.15; 0.83, 0.35]};
% 
% figure()
% hold on
% Gp = compositeGridPEBI([30,30], [1, 1], 'fracLines', l, 'padding', 1,...
%                        'fracGridSize', 0.05, 'circleFactor', 0.6);
% 
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight
% hold on
% plotFault(Gp)
% for i = 1:numel(l)
%   line = l{i};
%   %plot(line(:, 1), line(:, 2),'r');
% end



% %% Fault loosing points due to curvature.
% close all; clear
% 
% x = linspace(0.05,0.45);
% faultLine = {[(x+0.24)', 0.3*sin(2.0*pi*x)'+0.3]};
%             
% 
% Gp = compositeGridPEBI([15,15], [1, 1], 'faultLines', faultLine, 'padding', 1, ...
%                         'faultGridSize',0.09, 'mlqtMaxLevel', 2, ...
%                         'mlqtLevelSteps',[0.07,0.03]');
% 
% figure()
% hold on
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% %plotFault(Gp)
% %plotWells(Gp)
% for i = 1:numel(faultLine)
%   line = faultLine{i};
%   if size(line,1) == 1
%       plot(line(1,1), line(1,2),'.r', 'markersize', 8);
%   end
%   plot(line(:, 1), line(:, 2),'r');
% end

%%
% % Test of well grid points
% close all
% 
% x = linspace(0.2,0.8);
% wellLine = {[0.7105,0.1842], ...
%             [0.2,0.1842], ...
%             [x', 0.5*sin(pi*x)'+0.2]};%, ...
%             %[0.3,0.3;0.7,0.8]}
%             
% 
% Gp = compositeGridPEBI([1/19,0.02,1/40], [1, 1], 'wellLines', wellLine, 'padding', 1, ...
%                         'wellGridSize',0.02, 'mlqtMaxLevel', 2, ...
%                         'mlqtLevelSteps',[0.07,0.03]');
% 
% figure()
% hold on
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% %plotFault(Gp)
% %plotWells(Gp)
% for i = 1:numel(wellLine)
%   line = wellLine{i};
%   if size(line,1) == 1
%       plot(line(1,1), line(1,2),'.r', 'markersize', 8);
%   end
%   plot(line(:, 1), line(:, 2),'r');
% end
% 
% 
% 
%% Complex wells intersecting
% close all
% 
% x = linspace(0.2,0.8);
% wellLine = {[0.5,0.2; 0.5,0.3;0.47,0.4;0.4,0.5; 0.33,0.6;0.26,0.7], ...
%             [0.5,0.3;0.53,0.4;0.58,0.5],...
%             [0.5,0.45;0.5,0.55;0.45,0.65;0.4,0.75;0.38,0.85],...
%             [0.5,0.55;0.55,0.65;0.6,0.75;0.62,0.85]};
%                         
% 
% Gp = compositeGridPEBI(1/19, [1, 1], 'wellLines', wellLine,...
%                       'mlqtMaxLevel', 2, 'mlqtLevelSteps',[0.07,0.035]');
% 
% figure()
% hold on
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% %plotFault(Gp)
% %plotWells(Gp)
% for i = 1:numel(wellLine)
%   line = wellLine{i};
%   if size(line,1) == 1
%       plot(line(1,1), line(1,2),'.r', 'markersize', 8);
%   end
%   plot(line(:, 1), line(:, 2),'r');
% end

%% wells Intersecting fracture
% close all
% 
% wellLine = {[0.6,0.2;0.65,0.6],...        
%             [0.3,0.3;0.7,0.8],...
%             [0.6,0.2;0.85,0.4],...
%             [0.15,0.7;0.4,0.7]};
%         
% fracture = {[0.2,0.8;0.8,0.2]};
%       
% Gp = compositeGridPEBI([1/24,1/26/2,1/24/sqrt(2)], [1, 1], ...
%                        'wellLines', wellLine,'faultLines',fracture,...
%                         'circleFactor', 0.6);
% 
% figure()
% hold on
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% plotFault(Gp)
% plotWells(Gp)

%% Wells intersecting fracture multilevel
% close all
% 
% wellLine = {[0.6,0.2;0.65,0.6],...        
%             [0.3,0.3;0.7,0.8],...
%             [0.6,0.2;0.85,0.4],...
%             [0.15,0.7;0.4,0.7]};
%         
% fracture = {[0.2,0.8;0.8,0.2]};
%       
% Gp = compositeGridPEBI([1/24,0.006,0.03], [1, 1], 'wellLines', wellLine,'faultLines',fracture,...
%                          'circleFactor', 0.6,'padding', 1,...
%                          'mlqtMaxLevel', 2, ...
%                         'mlqtLevelSteps',[0.06,0.02]');
% 
% figure()
% hold on
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% plotFault(Gp)
% plotWells(Gp)

