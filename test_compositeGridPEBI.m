
% Multiple fractures. From Fung et.al 15.
close all
x = 0.2:0.05:0.8;
l = {[0.65,0.1;0.65,0.926],...
     [0.2,0.175; 0.875,0.875], ...
     [0.2,0.925; 0.9,0.125], ...
     [0.45,0.15; 0.83, 0.35]};

figure()
hold on
Gp = compositeGridPEBI([30,30], [1, 1], 'lines', l, 'padding', 1);

plotGrid(Gp, 'faceColor', 'none')
axis equal tight
hold on
plotFault(Gp)
%for i = 1:numel(l)
%   line = l{i};
%   plot(line(:, 1), line(:, 2),'r');
%end

%%
%Close up on one iregular fracture.
% x = [.114, 0.1432, 0.2341, 0.3477, 0.4795, 0.6205, 0.7705, 0.8695, 0.9014];
% y = [0.0991, 0.1050, 0.1574, 0.2391, 0.35048, 0.4006, 0.3956, 0.3786, 0.3515];
% l = {[x',y']};    
% 
% figure()
% hold on
% Gp = compositeGridPEBI([20,10], [1, 0.5], 'lines', l, 'padding', 1);
% 
% plotGrid(Gp, 'faceColor', 'none')
% axis equal tight
% hold on
% plotFault(Gp)