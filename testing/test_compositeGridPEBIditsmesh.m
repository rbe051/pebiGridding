% %% Single fracture
% close all; clear
% 
% wellLine = {}%{[0.2,0.2;0.67,0.2]}%,...        
% %            [0.3,0.3;0.7,0.8],...
% %            [0.6,0.2;0.85,0.4],...
% %            [0.15,0.7;0.4,0.7]};
%         
% fracture = {[0.2,0.8;0.8,0.2]};
%       
% Gp = compositeGridPEBIdistmesh([0.2,-1,-1], [1, 1], 'wellLines', wellLine,'faultLines',fracture,...
%                         'faultGridSize', 0.03, 'circleFactor', 0.6,'padding', 1,...
%                         'mlqtMaxLevel', 0, ...
%                         'mlqtLevelSteps',[0.06,0.02]');
% 
% figure()
% hold on
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% %plotFault(Gp)
% %plotWells(Gp)



% % % Multiple fractures. From Fung et.al 15.
% close all; clear all;
% x = 0.2:0.05:0.8;
% l = {[0.65,0.1;0.65,0.926],...is
%      [0.2,0.175; 0.875,0.875], ...
%      [0.2,0.925; 0.9,0.125], ...
%      [0.45,0.15; 0.83, 0.35]};
% 
% Gp = compositeGridPEBIdistmesh([1/30,-1,-1], [1, 1], 'faultLines', l, 'padding', 1);
% 
% figure()
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight
% hold on
% %plotFault(Gp)
% for i = 1:numel(l)sqrtsqrt
%   line = l{i};
%   plot(line(:, 1), line(:, 2),'r');
% end


% %%Complex faults intersecting.
% close all
% 
% x1 = linspace(0.35,0.2,10);
% y1 = [0.2,0.25,0.3,0.34,0.40,0.5,0.6,0.65,0.7,0.85];
% x2 = linspace(0.2,0.7,10);
% y2 = [0.45,0.4,0.4,0.38,0.35,0.35,0.4,0.5,0.6,0.73];
% x3 = linspace(0.1,0.9,10);
% y3 = [0.8,0.85,0.85,0.9,0.9,0.85,0.8,0.70,0.7,0.8];
% x4 = linspace(0.8,0.55,10);
% y4 = [0.2,0.25,0.3,0.34,0.40,0.5,0.6,0.65,0.7,0.85];
% x5 = linspace(0.2,0.7,10);
% y5 = [0.6,0.57,0.56,0.52,0.50,0.5,0.53,0.57,0.59,0.63];
% 
% 
% l = {[x2',y2'], [x3',y3'],  [x5',y5'],[x1',y1'],[x4',y4']};    
% 
% %l = {[0.4,0.2;0.8,0.8],[0.5,0.25;0.6,0.9],[0.4,0.7;0.8,0.65]};    
% 
% Gp = compositeGridPEBIdistmesh([1.155/50,-1,1/50], [1, 1], 'faultLines', l, 'padding', 0,...
%                         'fullFaultEdge', 1, 'circleFactor', 0.6);
% 
% figure()
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% %plotFault(Gp)
% 
% for i = 1:numel(l)
%   line = l{i};
%   plot(line(:, 1), line(:, 2),'r');
% end

%
% %Close up on one iregular fracture.
% close all; clear all
% x = [0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8];
% y = [0.25,0.26,0.28,0.34,0.39,0.4,0.39,0.36];
% l = {[x',y']};    
% 
% Gp = compositeGridPEBIdistmesh([0.1,-1,-1], [1, 0.75], 'faultLines', l, 'circleFactor', 0.6);
% 
% figure
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight off
% hold on
% plotFault(Gp)
% % %plotFault(Gp)



%%Test of well grid points. Sine cureve and two vertical wells
close all

x = linspace(0.2,0.8);
wellLine = {[0.7105,0.1842], ...
            [0.2,0.1842], ...
            [x', 0.5*sin(pi*x)'+0.2]};%, ...
            %[0.3,0.3;0.7,0.8]}
            

Gp = compositeGridPEBIdistmesh(1/20, [1, 1], 'wellLines', wellLine, ...
                               'wellGridFactor', 0.25);

figure()
hold on
plotGrid(Gp, 'faceColor', 'none')
axis equal tight off
hold on
%plotFault(Gp)

for i = 1:numel(wellLine)
  line = wellLine{i};
  if size(line,1) == 1
      plot(line(1,1), line(1,2),'.r', 'markersize', 8);
  end
  plot(line(:, 1), line(:, 2),'r');
end
plotWells(Gp)



