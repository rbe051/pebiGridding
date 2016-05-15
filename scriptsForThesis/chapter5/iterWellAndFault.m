clc; clear; close all


%% Create a well and fault Path

well  = {[0.2,0.2; 0.8,0.8], ...
         [0.4,0.4; 0.5,0.8]};
fault = {[0.2,0.7; 0.8,0.5]};

%% Plot
figure(1); hold on
col = get(gca,'ColorOrder');
plotLinePath(well,'color',col(1,:));
plotLinePath(fault,'color',col(2,:));
axis equal
axis([0,1,0,1])


%% Split well and fault at intersection
[sWell, ~, wfCut] = splitAtInt(well, fault);
[sFault,~, fwCut] = splitAtInt(fault, well);


%% Set gridSizes
df = 0.05;
dw = 0.05;
dr = 0.05;

F = createFaultGridPoints(sFault,df,'fwCut', fwCut);
fSites = F.f.pts;
fGs    = F.f.Gs;
[wSites, wGs] = createWellGridPoints(sWell, dw, 'wfCut', wfCut);

%% Plot
figure(2); hold on
col = get(gca,'ColorOrder');
plotLinePath(well,'color',col(1,:));
plotLinePath(fault,'color',col(2,:));
axis equal
axis([0,1,0,1])

plot(wSites(:,1), wSites(:,2),'.','color',col(1,:),'markersize',10)
plot(fSites(:,1), fSites(:,2),'.','color',col(2,:),'markersize',10)


%% Create reservoir sites.
[X,Y]  = meshgrid(0:dr:1);
rSites = [X(:),Y(:)];
rSites = removeConflictPoints2(rSites, [fSites;wSites], [fGs;wGs]);

%% Plot
figure(3); hold on
col = get(gca,'ColorOrder');
plotLinePath(well,'color',col(1,:));
plotLinePath(fault,'color',col(2,:));
axis equal
axis([0,1,0,1])

plot(wSites(:,1), wSites(:,2),'.','color',col(1,:),'markersize',10)
plot(fSites(:,1), fSites(:,2),'.','color',col(2,:),'markersize',10)

plot(rSites(:,1), rSites(:,2),'.k','markersize',10)


%% Create grid
Gt = triangleGrid([fSites;wSites;rSites]);
G = pebi(Gt);

%% Plot
figure(4)
plotGrid(G,'facecolor','none');
axis equal
axis([0,1,0,1])


%% Save plots
figure(1)
print('../../../../master/thesis/fig/ch05/interFaultWell1','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch05/interFaultWell2','-depsc')
figure(3)
print('../../../../master/thesis/fig/ch05/interFaultWell3','-depsc')
figure(4)
print('../../../../master/thesis/fig/ch05/interFaultWell4','-depsc')