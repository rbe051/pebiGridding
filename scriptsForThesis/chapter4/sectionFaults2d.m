clc; clear all; close all


%% Create a fault
fault    = [0.2,0.8; 0.5,0.5; 0.8,0.5];
dFault   = 0.1;
circFact = 0.6;

%% Plot Fault
figure()
hold on
plot(fault(:,1), fault(:,2));
axis equal tight
axis([0.1,0.9,0.4,0.9])
%print('../../../../master/thesis/fig/ch04/fault2DexampleImpJustFault','-depsc')
%% Interpolate fault
linesDist = sqrt(sum(diff(fault,[],1).^2,2));
linesDist = [0; linesDist]; % add the starting point
cumDist   = cumsum(linesDist);
dt        = cumDist(end)/ceil(cumDist(end)/dFault);
newPtsT   = 0:dt:cumDist(end);
CC        = interp1(cumDist, fault, newPtsT);

%% Plot interpolation:
figure()
hold on
yellow =  [0.9290    0.6940    0.1250];
plot(fault(:,1), fault(:,2));
plot(CC(:,1),CC(:,2),'.','color',yellow,'markersize',15);
axis equal tight
axis([0.1,0.9,0.4,0.9])
%print('../../../../master/thesis/fig/ch04/fault2DexampleImpInter','-depsc')

%% Calculate circle radiuses


d  = sqrt(sum((CC(2:end,:)-CC(1:end-1,:)).^2, 2));
CR = circFact*[d(1); (d(1:end-1) + d(2:end))/2; d(end)];
            
            
%% Calculate midpoints and circle offsets
a = (d.^2 - CR(2:end).^2 + CR(1:end-1).^2)./(2*d);
h = sqrt(CR(1:end-1).^2 - a.^2);
    
% Calculate normal vectors
n1 = (CC(2:end,:)-CC(1:end-1,:))./repmat(d,1,2); %Unit vector
n2 = [-n1(:, 2), n1(:,1)];                       %Unit normal


%% Calculate intersections
lPts = CC(1:end-1,:) + bsxfun(@times, a, n1)  ...
                     + bsxfun(@times, h,  n2);
rPts = CC(1:end-1,:) + bsxfun(@times, a, n1)  ...
                     - bsxfun(@times, h,  n2);
pts = [lPts; rPts];

%% Plotting
figure(); hold on
yellow =  [0.9290    0.6940    0.1250];

theta= linspace(0,2*pi)';
for i = 1:size(CC,1)
  X = repmat(CC(i,:),100,1) + repmat(CR(i),100,2).*[cos(theta), sin(theta)];
  plot(X(:,1),X(:,2),'k')
end
plot(fault(:,1), fault(:,2));
plot(pts(:,1),pts(:,2),'.','markersize',15);
plot(CC(:,1),CC(:,2),'.','color',yellow,'markersize',15);
axis equal tight
axis([0.1,0.9,0.4,0.9])
%print('../../../../master/thesis/fig/ch04/fault2DexampleImp','-depsc')

%% Cartesian background-grid
[X,Y]   = meshgrid(0:dFault:1);
backPts = [X(:),Y(:)];
CRrep   = repmat(CR',size(backPts,1),1);
removed = any(pdist2(backPts,CC)<CRrep,2);
backPts = backPts(~removed,:);

%% Plot all seeds
theta= linspace(0,2*pi)';
figure(); hold on
for i = 1:size(CC,1)
  X = repmat(CC(i,:),100,1) + repmat(CR(i),100,2).*[cos(theta), sin(theta)];
  plot(X(:,1),X(:,2),'k')
end
plot(fault(:,1), fault(:,2));
plot(pts(:,1),pts(:,2),'.','markersize',15);
plot(CC(:,1),CC(:,2),'.','color',yellow,'markersize',15);
plot(backPts(:,1), backPts(:,2),'.k','markersize',15);
axis equal tight
axis([0.1,0.9,0.4,0.9])
%print('../../../../master/thesis/fig/ch04/fault2DexampleImpAllSeeds','-depsc')

%% Create full grid
Gt = triangleGrid([pts;backPts]);
G  = pebi(Gt);
figure()
plotGrid(G,'facecolor','none')
axis equal tight
axis([0.1,0.9,0.4,0.9])
%print('../../../../master/thesis/fig/ch04/fault2DexampleImpGrid','-depsc')









