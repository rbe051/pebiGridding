clear; close all
%addpath('../../voronoi2D/')
%% set well and fault
f = {[0.2,0.4;0.8,0.6]};
w = {[0.5,0.2;0.5,0.8]};

%% Split at int
[fs,~,fwCut] = splitAtInt(f,w);
[ws,~,wfCut] = splitAtInt(w,f);

%% create points
df = 0.05;
dw = 0.05;
F = createFaultGridPoints(fs,df,'fwCut',fwCut);

bisectPnt = (df.^2 - (0.6*df).^2 + (0.6*df).^2)...
                ./(2*df);
faultOffset = sqrt((0.6*df).^2 - bisectPnt.^2);
sePtn = (1.0+faultOffset/dw)*[wfCut==2|wfCut==3, wfCut==1|wfCut==3];

[wSites,wGs] = createWellGridPoints(ws,dw,'sePtn',sePtn);

%% Create background grid
[X,Y]  = meshgrid(0:df*0.8:1);
rSites = [X(:),Y(:)];
rSites = removeConflictPoints2(rSites, [F.f.pts;wSites],[F.f.Gs;wGs]);

%%
pts = [F.f.pts;wSites;rSites];
Gt = triangleGrid(pts);
G  = pebi(Gt);
%% Plot
color = get(gca,'ColorOrder');
figure(1); hold on
theta = linspace(0,2*pi)';
for j = 1:size(F.c.CC,1)
X = repmat(F.c.CC(j,:),100,1) + repmat(F.c.R(j),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2));
end
axis equal
plotLinePath(f,'color',color(2,:))
plotLinePath(w,'color',color(1,:));
plot(F.f.pts(:,1), F.f.pts(:,2),'.','color',color(2,:),'markersize',15)
plot(wSites(:,1), wSites(:,2),'.','color',color(1,:),'markersize',15)
axis([0.35,0.65,0.35,0.65])
axis off

figure(2); hold on
plotGrid(G,'facecolor','none')
plotLinePath(f,'color',color(2,:))
plotLinePath(w,'color',color(1,:));
plot(F.f.pts(:,1), F.f.pts(:,2),'.','color',color(2,:),'markersize',15)
plot(wSites(:,1), wSites(:,2),'.','color',color(1,:),'markersize',15)
f = F.l.f([F.l.fPos([wfCut==2|wfCut==3;false]);F.l.fPos([false;wfCut==1|wfCut==3])-1]);
plot(F.f.pts(f,1), F.f.pts(f,2),'.','color',color(1,:),'markersize',15)
axis([0.35,0.65,0.35,0.65])
axis off


%% Save
figure(1);
print('../../../../master/thesis/fig/ch04/wellFaultIntersect','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch04/wellFaultIntersectGrid','-depsc')