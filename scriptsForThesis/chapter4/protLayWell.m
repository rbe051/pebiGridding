clear; close all

%% create well
w = {[0.3,0.3;0.7,0.7]};

%% Set protection distance function
gS = 0.2;
d = {@(x) 0.7*gS*ones(size(x,1),1)};

[wPts,wGs,protP,pGs] = createWellGridPoints(w,gS,'protLayer',true,'protD',d);
G = pebiGrid(gS,[1,1],'wellLines',w,'wellGridFactor',1,'protLayer',true,'protD',d);

%% Plot
close all
figure(); hold on
color = get(gca,'ColorOrder');
plotLinePath(w,'color',color(1,:))

plot(wPts(:,1), wPts(:,2),'.','color',color(1,:),'markersize',20)
plot(protP(:,1), protP(:,2),'.','color',color(4,:),'markersize',20)
axis equal off tight

figure(2); hold on
plotGrid(G,'facecolor','none')
plot(wPts(:,1), wPts(:,2),'.','color',color(1,:),'markersize',20)
plot(protP(:,1), protP(:,2),'.','color',color(4,:),'markersize',20)
axis equal off tight

%% Save
figure(1)
print('../../../../master/thesis/fig/ch04/inkscape/protLayer','-dsvg')
figure(2)
print('../../../../master/thesis/fig/ch04//inkscape/protLayerGrid','-dsvg')