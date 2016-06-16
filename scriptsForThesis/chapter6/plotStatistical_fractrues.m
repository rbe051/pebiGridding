clear; 

%% load
load('statistical_fractures_edgeCentered.mat')
% Rotate
l = cellfun(@(c) [c(:,2), c(:,1)], l,'un',false);
G.nodes.coords = [G.nodes.coords(:,2), G.nodes.coords(:,1)];
%% Plot grid
pdims = max(G.nodes.coords);
bdr = [0,0;0,pdims(2);pdims;pdims(1),0;0,0];
ax = [91,98,22,28];
box1 = [ax(1),ax(3);ax(1),ax(4); ax(2),ax(4); ax(2), ax(3); ax(1),ax(3)];
ax2 = [56,62,21,27];
box2 = [ax2(1),ax2(3);ax2(1),ax2(4); ax2(2),ax2(4); ax2(2), ax2(3); ax2(1),ax2(3)];



color = get(gca,'ColorOrder');
close all
figure(1)
plotLinePath(l,'color','r','linewidth',1.5);
plot(bdr(:,1), bdr(:,2),'k--')
axis equal tight off
figure(2); hold on
plotGrid(G,'facecolor','none');
plot(box1(:,1), box1(:,2),'--','color',color(4,:),'linewidth',2);
plot(box2(:,1), box2(:,2),'--','color',color(5,:),'linewidth',2)
axis equal tight off

figure(3); hold on
plotGrid(G,'facecolor','none')
plotLinePath(l,'r--','linewidth',1)
plot(box1(:,1), box1(:,2),'--','color',color(4,:),'linewidth',2)
axis tight off
axis(ax)

figure(4); hold on

plotGrid(G,'facecolor','none')
plotLinePath(l,'r--','linewidth',1)
plot(box2(:,1), box2(:,2),'--','color',color(5,:),'linewidth',2)
axis tight off
axis(ax2)
%% Save
figure(1)
pause(0.1)
print('../../../../master/thesis/fig/ch06/statisticalFracturesOnlyFrac','-depsc','-painters')


figure(2)
pause(0.1)
print('../../../../master/thesis/fig/ch06/statisticalFracturesGrid','-depsc','-painters')

figure(3)
pause(0.1)
print('../../../../master/thesis/fig/ch06/statisticalFracturesZoom1','-depsc','-painters')

figure(4)
pause(0.1)
print('../../../../master/thesis/fig/ch06/statisticalFracturesZoom2','-depsc','-painters')