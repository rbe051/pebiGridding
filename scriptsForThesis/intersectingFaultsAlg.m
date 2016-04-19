clc; clear; close all; addpath ../voronoi2D/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fault intersecting fault
fl = {[0.5,0.5;0.8,0.8],[0.5,0.5;0.8,0.5]};

figure()
hold on
for i = 1:numel(fl)
  f = fl{i};
  plot(f(:,1), f(:,2),'k');
end

G  = compositePebiGrid(1/10,[1,1],'faultLines', fl);
% Paste the following into createFaultGridPoitns., -> fixIntersection(F)
%% Before merge
   blue  = [0    0.4470    0.7410];
   theta = linspace(0,2*pi)';
   cc = circ(:);
  for i = 1:size(F.c.R,1)
    X = repmat(F.c.CC(i,:),100,1) + repmat(F.c.R(i),100,2).*[cos(theta), sin(theta)];
    plot(X(:,1), X(:,2),'color',blue)
  end
  for i = 1:size(cc)
      plot(F.f.pts(F.c.f(F.c.fPos(cc(i)):F.c.fPos(cc(i)+1)-1),1),F.f.pts(F.c.f(F.c.fPos(cc(i)):F.c.fPos(cc(i)+1)-1),2),'r.','markersize',20)
  end
 axis([0.46,0.57,0.46,0.57]);
  axis off
%% after merge
  figure();
  hold on
     c = circ(:);
  for i = 1:size(F.c.R)
    X = repmat(F.c.CC(i,:),100,1) + repmat(F.c.R(i),100,2).*[cos(theta), sin(theta)];
    plot(X(:,1), X(:,2),'color',blue)
  end
  for i = 1:size(cc)
      plot(F.f.pts(F.c.f(F.c.fPos(cc(i)):F.c.fPos(cc(i)+1)-1),1),F.f.pts(F.c.f(F.c.fPos(cc(i)):F.c.fPos(cc(i)+1)-1),2),'r.','markersize',20)
  end

  axis([0.46,0.57,0.46,0.57]);
  axis off
  %%

for i = 1:numel(fl)
  f = fl{i};
  plot(f(:,1), f(:,2),'k');
end

figure()

plotGrid(G, 'FaceColor', 'none')

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
  plot(w(:,1), w(:,2),'color',blue)
end

G  = compositePebiGrid(1/10,[1,1],'faultLines', fl,'wellLines',wl);
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

