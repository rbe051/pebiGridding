%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fault intersecting well
clc; clear; close all
fl = {[0.2,0.4;0.8,0.6]};
wl = {[0.5,0.2;0.5,0.8]};

figure()
hold on
blue = [0    0.4470    0.7410];
red  = [0.8500  0.3250  0.0980];
for i = 1:numel(fl)
  f = fl{i};
  w = wl{i};
  plot(f(:,1), f(:,2),'color',red);
  %plot(w(:,1), w(:,2),'color',blue)
end

G  = compositePebiGrid(1/10,[1,1],'faultLines', fl,'wellLines',wl);
plotGrid(G,'facecolor','none')
%% PASTE INTO END OF COMPOSITEPEBIGRID()
blue = [0    0.4470    0.7410];
red  = [0.8500  0.3250  0.0980];
theta = linspace(0,2*pi)';
for i = 1:size(F.c.CC)
X = repmat(F.c.CC(i,:),100,1) + repmat(F.c.R(i),100,2).*[cos(theta), sin(theta)];
plot(X(:,1), X(:,2),'k')
end
plot(F.f.pts(:,1), F.f.pts(:,2),'.','color',red,'markersize',20);
plot(wellPts(:,1), wellPts(:,2),'.','color',blue,'markersize',20),
axis([0.44,0.56,0.44,0.56])
axis off

