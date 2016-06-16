clear; close all
  

%% create fault
fault = {[0.45,0.2;0.55,0.8],...
         [0.55,0.2;0.45,0.8]};
%fault{1} = [fault{1}(:,2),fault{1}(:,1)];
%fault{2} = [fault{2}(:,2),fault{2}(:,1)];
fGs = 1/12;
[sFault,fCut,~] = splitAtInt(fault,{});

F = createFaultGridPoints(sFault,fGs,'fCut',fCut);

color = get(gca,'ColorOrder');
[~,I] = max(diff(F.c.lPos));
theta = linspace(0,2*pi)';
X = repmat(F.c.CC(I,:),100,1) + repmat(F.c.R(I),100,2).*[cos(theta), sin(theta)];
plot(X(:,1), X(:,2),'color',color(5,:))
plot(F.f.pts(:,1), F.f.pts(:,2),'.','color',color(2,:),'markersize',15);

%% Create background grid
[X,Y]  = meshgrid(0:fGs*0.8:1);
rSites = [X(:),Y(:)];
rSites = removeConflictPoints2(rSites, F.f.pts,F.f.Gs);

%%
pts = [F.f.pts;rSites];
Gt = triangleGrid(pts);
G  = pebi(Gt);

%%
figure(); hold on
plotGrid(G,'facecolor','none')
plotLinePath(fault,'--','color',color(2,:));

%% Save 
for i = 1:3
  figure(i)
  axis equal off;
  axis([0.35,0.65,0.3,0.7])
  fPath = strcat('../../../../master/thesis/fig/ch04/inkscape/mergeCircle',num2str(i));
  print(fPath,'-dsvg')
end
  
