clear; close all



fault = {[0,0; 0.5,0.0; 0.5+0.5/sqrt(2),0.5/sqrt(2)]};
df = 1/2;

%% Circle factor 0.9
F = createFaultGridPoints(fault, df,'circleFactor',0.9);

color = get(gca,'colorOrder');
figure(); hold on
plotLinePath(fault,'color',color(2,:))

theta = linspace(0,2*pi)';
for j = 1:size(F.c.CC,1)
X = repmat(F.c.CC(j,:),100,1) + repmat(F.c.R(j),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'color','k');
end

plot(F.f.pts(:,1), F.f.pts(:,2),'.','color',color(2,:), 'markersize',20)
plot(F.c.CC(:,1),F.c.CC(:,2),'.k','MarkerSize',20)
axis equal off

%% Circle factor 0.6
F = createFaultGridPoints(fault, df,'circleFactor',0.6);
limx = get(gca,'xlim');
limy = get(gca,'ylim');

color = get(gca,'colorOrder');
figure(); hold on
plotLinePath(fault,'color',color(2,:))

theta = linspace(0,2*pi)';
for j = 1:size(F.c.CC,1)
X = repmat(F.c.CC(j,:),100,1) + repmat(F.c.R(j),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2),'color','k');
end

plot(F.f.pts(:,1), F.f.pts(:,2),'.','color',color(2,:), 'markersize',20)
plot(F.c.CC(:,1),F.c.CC(:,2),'.k','MarkerSize',20)
axis equal off
axis([limx,limy])