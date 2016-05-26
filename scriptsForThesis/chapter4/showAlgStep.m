
addpath ../../voronoi2D/

%% Create a fault and well
clc; clear all; close all
fault    = {[0.2,0.8; 0.8,0.2]};
well     = {[0.3,0.2;0.3,0.35;0.35,0.5;0.45,0.6;0.5,0.65;0.8,0.8],...
            [0.5,0.8;0.8,0.5]};
dFault   = 0.1;
circFact = 0.6;

offset = 0.05;
figure(); hold on
orange = [0.8500,    0.3250,    0.0980];
blue   = [0 ,   0.4470,    0.7410];
axis off equal tight
%axis([-offset,offset,-offset,offset])
for i = 1:numel(well)
  plot(well{i}(:,1), well{i}(:,2),'color',blue)
end
for i = 1:numel(fault)
  plot(fault{i}(:,1), fault{i}(:,2),'color',orange)
end


boxX = [-offset,-offset,1+offset,1+offset,-offset];
boxY = [-offset,1 + offset,1+offset,-offset,-offset];
plot(boxX, boxY, '--','color',[0.25,0.25,0.25])


G = compositePebiGrid(dFault,[1,1],'wellLines',well,'faultLines',fault,...
                      'wellGridFactor',0.5,'faultGridFactor',1/sqrt(2), ...
                      'mlqtMaxLevel', 1);

figure()                    
plotGrid(G,'facecolor','none')
axis equal tight off
print('../../../../master/thesis/fig/ch04/showAlgoGrid','-depsc')

%% Plotting: paste before remove conflict points
color = get(gca,'colororder');
orange = [0.8500,    0.3250,    0.0980];
blue   = [0 ,   0.4470,    0.7410];
plot(F.f.pts(:,1),F.f.pts(:,2),'.','color',orange,'markersize',15);
plot(wellPts(:,1),wellPts(:,2),'.','color',blue,'markersize',15);
f = F.l.f([F.l.fPos([wfCut==2|wfCut==3;false]);F.l.fPos([false;wfCut==1|wfCut==3])-1]);
plot(F.f.pts(f,1), F.f.pts(f,2),'.','color',color(1,:),'markersize',15)
print('../../../../master/thesis/fig/ch04/showAlgoWellAndFaultPts','-depsc')
a = plot(resPtsInit(:,1), resPtsInit(:,2),'.k','markersize',15);
print('../../../../master/thesis/fig/ch04/showAlgoWellAndFaultPtsAndResPts','-depsc')
delete(a)
a = plot(resPts(:,1), resPts(:,2),'.k','markersize',15);
print('../../../../master/thesis/fig/ch04/showAlgoWellAndFaultPtsAndResPtsRef','-depsc')


%% Plotting: paste after remove conflict points
delete(a)
a = plot(resPts(:,1), resPts(:,2),'.k','markersize',15);
print('../../../../master/thesis/fig/ch04/showAlgoWellAndFaultPtsAndResPtsRem','-depsc')