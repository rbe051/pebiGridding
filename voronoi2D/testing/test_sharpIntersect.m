%% Faults intersecting

%close all

l = {[0.2,0.48;0.8,0.52], ...
     [0.2,0.52;0.8,0.48], ...
     [0.5,0.2; 0.5,0.8]};

%l = {[0.4,0.2;0.8,0.8],[0.5,0.25;0.6,0.9],[0.4,0.7;0.8,0.65]};     

G = compositePebiGrid(1.155/50, [1, 1], 'faultLines', l, 'faultGridFactor', 1/1.155,...
                        'fullFaultEdge', 1, 'circleFactor', 0.6);

figure()
plotGrid(G,'facecolor','none')
axis equal tight off
hold on


