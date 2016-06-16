addpath ../../voronoi2D/

     
% Paste the following into createFaultGridPoitns., -> fixIntersection(F)
%% paste into createFaultGridPoints before if statement about merging circles
color = get(gca,'ColorOrder');
theta = linspace(0,2*pi)';
for j = 1:size(F.c.CC,1)
X = repmat(F.c.CC(j,:),100,1) + repmat(F.c.R(j),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'k');
end
for j = 1:size(circ,1)
X = repmat(F.c.CC(circ(j,3),:),100,1) + repmat(F.c.R(circ(j,3)),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'color',color(1,:));
X = repmat(F.c.CC(circ(j,1),:),100,1) + repmat(F.c.R(circ(j,1)),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'color','r');
X = repmat(F.c.CC(circ(j,2),:),100,1) + repmat(F.c.R(circ(j,2)),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'color','r');
end
plot(F.f.pts(:,1), F.f.pts(:,2),'.','color',color(2,:),'markersize',25)
plot(F.c.CC(:,1), F.c.CC(:,2),'.k','markersize',15)
axis off equal
axis([0.4,0.7,0.4,0.7])
%%

 %% after merge
  figure();
  hold on
for j = 1:size(F.c.CC,1)
X = repmat(F.c.CC(j,:),100,1) + repmat(F.c.R(j),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'k');
end
for j = 1:size(circ,1)
X = repmat(F.c.CC(circ(j,3),:),100,1) + repmat(F.c.R(circ(j,3)),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'color',color(1,:));
X = repmat(F.c.CC(circ(j,1),:),100,1) + repmat(F.c.R(circ(j,1)),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'color','r');
X = repmat(F.c.CC(circ(j,2),:),100,1) + repmat(F.c.R(circ(j,2)),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'color','r');
end
plot(F.c.CC(:,1), F.c.CC(:,2),'.k','markersize',15)
plot(F.f.pts(:,1),F.f.pts(:,2),'.','color',color(2,:),'markersize',25)
axis off equal
axis([0.4,0.7,0.4,0.7])
  figure(3)
  plot(F.f.pts(:,1),F.f.pts(:,2),'.','color',color(2,:),'markersize',25)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fault intersecting fault
clc; clear; close all; 
fl = {[0.5,0.5;0.8,0.8],[0.5,0.5;0.8,0.5]};
color = get(gca,'ColorOrder')

hold on
for i = 1:numel(fl)
  f = fl{i};
  plot(f(:,1), f(:,2),'color',color(2,:));
end

G  = compositePebiGrid([1/10,1/10],[1,1],'faultLines', fl);

figure(2)

plotLinePath(fl,'color',color(2,:));


figure(3);
plotGrid(G, 'FaceColor', 'none')
axis([0.46,0.57,0.46,0.57]);
axis off equal
axis([0.4,0.7,0.4,0.7])



%% save
figure(1)
pause(0.1)
print('~/masterThesis/master/thesis/fig/ch04/intersectingFaultsMerge1','-depsc','-painters')
figure(2)
pause(0.1)
print('~/masterThesis/master/thesis/fig/ch04/intersectingFaultsMerge2','-depsc','-painters')
figure(3)
pause(0.1)
print('~/masterThesis/master/thesis/fig/ch04/intersectingFaultsMegeGrid','-depsc','-painters')
