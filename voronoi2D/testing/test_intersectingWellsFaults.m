clear; %close all

%% Faults intersecting

%close all

x5 = linspace(0.35,0.2,10);
y5 = [0.2,0.25,0.3,0.34,0.40,0.5,0.6,0.65,0.7,0.85];
x4 = linspace(0.2,0.7,10);
x4(8) = 0.56;
x4(9) = 0.65;
y4 = [0.45,0.4,0.4,0.38,0.35,0.35,0.4,0.48,0.615,0.73];
x3 = linspace(0.1,0.9,10);
y3 = [0.8,0.85,0.85,0.9,0.9,0.85,0.8,0.70,0.7,0.8];
x2 = linspace(0.8,0.55,8);
x2(7) = 0.62;
x2(6) = 0.62;
y2 = [0.2,0.25,0.3,0.34,0.40,0.5,0.7,0.85];
x1 = linspace(0.2,0.7,10);
y1 = [0.6,0.57,0.56,0.52,0.50,0.5,0.53,0.57,0.57,0.63];


l = {[x2',y2'], [x3',y3'],  [x5',y5'],[x1',y1'],[x4',y4']};    

%l = {[0.4,0.2;0.8,0.8],[0.5,0.25;0.6,0.9],[0.4,0.7;0.8,0.65]};     

G = compositePebiGrid(1.155/50, [1, 1], 'faultLines', l, 'faultGridFactor', 1/1.155,...
                        'fullFaultEdge', 1, 'circleFactor', 0.6);
Gp = pebiGrid(1.155/50, [1, 1], 'faultLines', l, 'faultGridFactor', 1/1.155,...
                        'circleFactor', 0.6);
figure()
plotGrid(G)                  
figure()
plotGrid(Gp, 'faceColor', 'none')
axis equal tight off
hold on


for i = 1:numel(l)
  line = l{i};
  plot(line(:, 1), line(:, 2),'r');
end


%% Wells Intersecting

x5 = linspace(0.35,0.2,10);
y5 = [0.2,0.25,0.3,0.34,0.40,0.5,0.6,0.65,0.7,0.85];
x4 = linspace(0.2,0.7,10);
x4(8) = 0.56;
x4(9) = 0.65;
y4 = [0.45,0.4,0.4,0.38,0.35,0.35,0.4,0.48,0.615,0.73];
x3 = linspace(0.1,0.9,10);
y3 = [0.8,0.85,0.85,0.9,0.9,0.85,0.8,0.70,0.7,0.8];
x2 = linspace(0.8,0.55,8);
x2(7) = 0.62;
x2(6) = 0.62;
y2 = [0.2,0.25,0.3,0.34,0.40,0.5,0.7,0.85];
x1 = linspace(0.2,0.7,10);
y1 = [0.6,0.57,0.56,0.52,0.50,0.5,0.53,0.57,0.57,0.63];


l = {[x2',y2'], [x3',y3'],  [x5',y5'],[x1',y1'],[x4',y4']};    

%l = {[0.4,0.2;0.8,0.8],[0.5,0.25;0.6,0.9],[0.4,0.7;0.8,0.65]};     

Gp = compositePebiGrid(1.155/20, [1, 1], 'wellLines', l, 'wellGridFactor', 1/1.155);
G = pebiGrid(1.155/20, [1, 1], 'wellLines', l, 'wellGridFactor', 1/1.155);

figure()
plotGrid(Gp)

figure()
plotGrid(G, 'faceColor', 'none')
axis equal tight off
hold on
%plotFault(Gp)

for i = 1:numel(l)
  line = l{i};
  plot(line(:, 1), line(:, 2),'r');
end
G = computeGeometry(G);
plot(G.cells.centroids(G.cells.tag,1),G.cells.centroids(G.cells.tag,2),'.','markersize',15)

%% Faults and wells intersecting
x1 = [0.1,.9];
y1 = [0.1,.9];
x2 = [0.1,.9];
y2 = [0.3,0.2];
lw = {[x1',y1'],[x2',y2']};

x3 = [0.5,0.5];
y3 = [0.2,0.8];
x4 = [0.2,0.8];
y4 = [0.8,0.6];
lf ={[x3',y3'], [x4',y4']};
  

Gp = compositePebiGrid(1.155/20, [1, 1], 'wellLines', lw, 'wellGridFactor', 1/4,...
                      'mlqtMaxLevel',2,...
                      'faultLines', lf,'faultGridFactor', 1/4);
% G  = pebiGrid(1.155/20,[1,1], 'wellLines', lw, 'wellGridFactor', 1/4,...
%              'wellRefinement',true,'epsilon',1/10,...
%                        'faultLines', lf,'faultGridFactor', 1/1.155);
G = computeGeometry(G);
Gp =  computeGeometry(Gp);
figure()
plotGrid(Gp, 'faceColor', 'none')
axis equal tight off
hold on
plot(Gp.cells.centroids(Gp.cells.tag,1), Gp.cells.centroids(Gp.cells.tag,2),'.','markersize',10)
plot(Gp.faces.centroids(Gp.faces.tag,1), Gp.faces.centroids(Gp.faces.tag,2),'.','markersize',10)
figure()
plotGrid(G)
axis equal tight off
hold on
plot(G.faces.centroids(G.faces.tag,1), G.faces.centroids(G.faces.tag,2),'.','markersize',10)
plot(G.cells.centroids(G.cells.tag,1), G.cells.centroids(G.cells.tag,2),'.','markersize',10)
