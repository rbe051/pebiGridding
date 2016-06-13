clear; close all
%%

load('showInnerOuterCellsElephant.mat');


%% Plot grids
color = rand(G.cells.num,3);
figure(1); hold on
axis equal off tight
figure(2); hold on
axis equal off tight
for i = 1:G.cells.num
  figure(1)
  plotGrid(G,i,'facecolor',color(i,:));
  figure(2)  
  plotGrid(Gc,i,'facecolor',color(i,:));
end
G = computeGeometry(G);
Gc = computeGeometry(Gc);
%% Plot sites
figure(1)
ptsInit = [ptsx,ptsy];
plot(pts(:,1), pts(:,2),'.k','markersize',10)
%plot(G.cells.centroids(:,1), G.cells.centroids(:,2),'ko')
figure(2)
plot(optPts(:,1), optPts(:,2),'.k','markersize',10)
%plot(Gc.cells.centroids(:,1), Gc.cells.centroids(:,2),'ko')


%% Save figure
figure(1)
print('../../../../master/thesis/fig/ch03/elephantInitial','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch03/elephantConverged','-depsc')