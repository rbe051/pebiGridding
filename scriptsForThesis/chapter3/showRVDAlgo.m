clear; close all
addpath ../
addpath ../../voronoi3D/
addpath ../../../vem/mat/
%% Creat domain
k = 10;
r = 1;
theta = linspace(0,2*pi*(k-1)/k,k)';
bound = [r*cos(theta), r*sin(theta)];
boundAtInf = bound*200;
bd = delaunayTriangulation(bound);


%% Create sites
n = 15;
rng(11);
pts = 1 - 2*rand(n,2);
pts = pts(~isnan(bd.pointLocation(pts)),:);
ptsAtInf= [-100,-100;
           -100, 100;
            100, 100;
            100,-100];
          
ptsAtInf = [pts;ptsAtInf];


%% Create RVD grid
G = clippedPebi2D(pts,bound);
G = computeGeometry(G);

%% Create Voronoi Grid
Ginf = clippedPebi2D(ptsAtInf, boundAtInf);
Ginf = computeGeometry(Ginf);          

%% Plot
figure(); hold on
plotGrid(Ginf,'facecolor','none')
plot([bound(:,1);bound(1,1)], [bound(:,2);bound(1,2)],'--k')
axis([-1.2,1.2,-1.2,1.2])
axis equal off
print('../../../../master/thesis/fig/ch03/notRVDshowAlgo','-depsc')

%%
figure()
plotGrid(G,'facecolor','none')
axis([-1.2,1.2,-1.2,1.2])
axis off equal

%%
R = G;
R.nodes.coords(:,1:2) = 1.00001*R.nodes.coords(:,1:2);

R = computeGeometry(R);
rFace = isnan(bd.pointLocation(R.faces.centroids));
rFace = find(rFace);

colorOrder = get(gca, 'ColorOrder');
nc = size(colorOrder,1);
cell = reshape(G.faces.neighbors(rFace,:)',[],1);
cell = cell(cell~=0);
[~,~,IC] = unique(cell);
color = colorOrder(mod(IC,nc-1)+1,:);

for i = 1:numel(rFace)
  c = color(i,:);
  plotFaces(G, rFace(i), 'facecolor','none','edgecolor',c,'linewidth',1)
end

print('../../../../master/thesis/fig/ch03/RVDshowAlgo','-depsc')