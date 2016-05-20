clear all; close all

%% Set fault
l = {[0.1,0.42; 0.4,.55; 0.7,0.65], ...
     [0.8,0.13; 0.6,0.4; 0.55,0.6],...
     [0.42,1.08; 0.45,0.9; 0.5,0.8; 0.58,0.6]};
   
gs = 1/30;
   
%% Create grid
G = compositePebiGrid(1/30,[1,1.15],'faultLines',l);

%% Plot grid
figure(1); hold on
plotGrid(G,'facecolor','none');


plotLinePath(l,'-','color','k','linewidth',1)

%% Plot box for zoom
color = get(gca,'ColorOrder');
boxX = [0.4,0.4,0.7,0.7,0.4];
boxY = [0.45,0.75,0.75,0.45,0.45];
plot(boxX, boxY, '--','color',color(3,:),'linewidth',2)%[0.25,0.25,0.25])
axis equal off tight

%% Plot zoomed
figure(2); hold on
plotGrid(G,'facecolor','none');
plotLinePath(l,'-','color','k','linewidth',1.5)
axis equal off tight
axis([0.4,0.7,0.4,0.7 + (0.7-0.4)*.15])

%% save
figure(1)
print('../../../../master/thesis/fig/ch06/exFromFung','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch06/exFromFungZoom','-depsc')