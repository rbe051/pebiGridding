close; clear all
%addpath ../../voronoi2D/

%%
color = get(gca,'ColorOrder');
%% create faults and wells
f = {[0.2,0.2;0.5,0.5;0.8,0.6], ...
     [0.2,0.8;0.5,0.5]};
w = {[0.8,0.2;0.7,0.4;0.3,0.8], ...
     [0.7,0.4;0.5,0.35]};


%% Set grid parameters
pdims = [1,1];

fGf  = 1/2;
wGf = 1/2;
rGs = 0.08;

mlqtLevels = 1;
eps = 1/5;

%% Create Grid
G = compositePebiGrid(rGs,pdims,'faultLines',f,         ...
                                'wellLines', w,         ...
                                'faultGridFactor', fGf, ...
                                'wellGridFactor',  wGf, ...
                                'mlqtMaxLevel',1);
Gp = pebiGrid(rGs,pdims,'faultLines',f,         ...
                        'wellLines', w,         ...
                        'faultGridFactor', fGf, ...
                        'wellGridFactor',  0.8*wGf, ...
                        'wellRefinement',true,...
                        'epsilon', eps);
                       
%%
G = computeGeometry(G);
Gp = computeGeometry(Gp);
%% plot
% Well and fault paths
figure(1);
plotLinePath(f,'color',color(2,:))
plotLinePath(w,'color',color(1,:))
axis equal
axis([0,1,0,1])
% grid 1
figure(2); hold on
plotGrid(G,'facecolor','none')
c = G.cells.tag;
fa = G.faces.tag;
plot(G.cells.centroids(c,1),G.cells.centroids(c,2),'.','color',color(1,:))
plotFaces(G,fa,'edgecolor',color(2,:))
axis equal off tight
% grid2
figure(3); hold on
plotGrid(Gp,'facecolor','none')
c = Gp.cells.tag;
fa = Gp.faces.tag;
plot(Gp.cells.centroids(c,1),Gp.cells.centroids(c,2),'.','color',color(1,:))
plotFaces(Gp,fa,'edgecolor',color(2,:))

axis equal off tight

%% Save
figure(1)
print('../../../../master/thesis/fig/ch05/showWrappFunc1','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch05/showWrappFunc2','-depsc')
figure(3)
print('../../../../master/thesis/fig/ch05/showWrappFunc3','-depsc')