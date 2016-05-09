clear; close all


%% Creat domain
k = 10;
r = 1;
theta = linspace(0,2*pi*(k-1)/k,k)';
bound = [r*cos(theta), r*sin(theta)];
bound = [bound, -0.01*ones(k,1); bound, 0.1*ones(k,1)];
boundAtInf = [bound(:,1:2)*200,bound(:,3)];
bound = delaunayTriangulation(bound);
boundAtInf = delaunayTriangulation(boundAtInf);

%% Create sites
n = 15;

pts = 1 - 2*rand(n,3);
pts(:,3) = 0.01*pts(:,3);
pts = pts(~isnan(bound.pointLocation(pts)),:);
ptsAtInf= [-100,-100,0;
           -100, 100,0;
            100, 100,0;
            100,-100,0];
          
ptsAtInf = [pts;ptsAtInf];



G = restrictedVoronoiDiagram(pts,bound);
G = computeGeometry(G);
Ginf = restrictedVoronoiDiagram(ptsAtInf, boundAtInf);
Ginf = computeGeometry(Ginf);          

ninf = bsxfun(@times, Ginf.faces.normals, (1-2*(Ginf.faces.neighbors(:,2)~=0)));
finf = sum(ninf(:,1:2).^2,2)<1e-10 & ninf(:,3)>0;
ng = bsxfun(@times, G.faces.normals, (1-2*(G.faces.neighbors(:,2)~=0)));
fg = sum(ng(:,1:2).^2,2)<1e-10 & ng(:,3)>0;

figure()
plotFaces(Ginf,finf,'facecolor','none')
axis([-1.2,1.2,-1.2,1.2])
axis equal off
view(0,-90)
print('../../../../master/thesis/fig/ch03/notRVDshowAlgo','-depsc')


figure()
plotFaces(G,fg,'facecolor','none')
axis([-1.2,1.2,-1.2,1.2])
axis off equal tight


R = G;
R.nodes.coords(:,1:2) = 1.00001*R.nodes.coords(:,1:2);
R.nodes.coords(:,3) = 0.5*R.nodes.coords(:,3);

R = computeGeometry(R);
rFace = isnan(bound.pointLocation(R.faces.centroids));
rFace = find(rFace);

colorOrder = get(gca, 'ColorOrder');
nc = size(colorOrder,1);
cell = reshape(G.faces.neighbors(rFace,:)',[],1);
cell = cell(cell~=0);
color = colorOrder(mod(cell,nc-1)+1,:);
G.nodes.coords(:,3) = 1.01*G.nodes.coords(:,3);
for i = 1:numel(rFace)
  c = color(i,:);
  plotFaces(G, rFace(i),'facecolor','none','edgecolor',c,'linewidth',1)
  
end
view(0,-90)

print('../../../../master/thesis/fig/ch03/RVDshowAlgo','-depsc')