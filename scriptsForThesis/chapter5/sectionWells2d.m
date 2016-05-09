clc; clear all; close all


%% Create a Well
dWell = 0.1;
n     = 20;
dT    = 2*pi/20;
theta = 2*pi*(0:1/n:1 -1/n)';
well  = repmat([0.5,0.5],n,1) ...
          + 0.3*[cos(theta), sin(theta)];

%% Plot Fault
figure()
hold on
plot([well(:,1);well(1,1)], [well(:,2);well(1,2)]);
axis equal tight
axis([0.1,0.9,0.1,0.9])
print('../../../../master/thesis/fig/ch04/well2DexampleImpJustWell','-depsc')
%% Interpolate fault
linesDist = sqrt(sum(diff(well,[],1).^2,2));
linesDist = [0; linesDist]; 
cumDist   = cumsum(linesDist);
dt        = cumDist(end)/...
              ceil(cumDist(end)/dWell);
newPtsT   = 0:dt:cumDist(end);
wPts   = interp1(cumDist, well, newPtsT);

%% Plot interpolation:
figure()
hold on
yellow = [0.9290    0.6940    0.1250];
blue   = [0    0.4470    0.7410];
plot([well(:,1);well(1,1)], [well(:,2);well(1,2)]);
plot(wPts(:,1),wPts(:,2),'.','color',blue,'markersize',15);
axis equal tight
axis([0.1,0.9,0.1,0.9])
print('../../../../master/thesis/fig/ch04/well2DexampleImpInter','-depsc')

%% Cartesian background-grid
[X,Y]   = meshgrid(0:dWell:1);
resPts = [X(:),Y(:)];
removed = any(pdist2(resPts,wPts)<dWell,2);
resPts = resPts(~removed,:);

%% Plot all seeds
figure(); hold on
plot([well(:,1);well(1,1)], [well(:,2);well(1,2)]);
plot(wPts(:,1),wPts(:,2),'.','color',blue,'markersize',15);
plot(resPts(:,1), resPts(:,2),'.k','markersize',15);
axis equal tight
axis([0.1,0.9,0.1,0.9])
print('../../../../master/thesis/fig/ch04/well2DexampleImpAllSeeds','-depsc')

%% Create full grid
Gt = triangleGrid([wPts;resPts]);
G  = pebi(Gt);
figure()
plotGrid(G,'facecolor','none')
axis equal tight
axis([0.1,0.9,0.1,0.9])
print('../../../../master/thesis/fig/ch04/well2DexampleImpGrid','-depsc')









