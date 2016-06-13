clear; close all
addpath ../
addpath ../../voronoi3D/
addpath ../../../vem/mat/
addpath ../../../vem/mat/VEM2D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Packman
%% Create domain
phi = pi/4;
k   = 30;
theta = linspace(phi, 2*pi - phi,k);
bdr   = [cos(theta)', sin(theta)';...
         0          ,          0];



%% Create seeds
n = 200;
shift = (2*pi-1.1*2*phi)*rand(n,1) + 1.1*phi;
r     = 0.9*sqrt(rand(n,1))+0.05;
pts = bsxfun(@times, [cos(shift),sin(shift)], r);
%% Create grid

G = clippedPebi2D(pts, bdr);
Gc = createCVD(pts,bdr);
%% plot

plotGrid(G)
axis equal off tight
figure()
plotGrid(Gc)
axis equal off tight
%% Save grids
save('showInnerOuterCellsPackMan','G','Gc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sheep
%% Load domain

mat = load('sheep.mat');
bdr = mat.P;
bdr = bdr/norm(bdr,'inf');

plot(bdr(:,1), bdr(:,2));
hold on
%% Create seeds.
n = 300;
x = linspace(min(bdr(:,1))+1/n,max(bdr(:,1))-1/n,round(sqrt(n)));
y = linspace(min(bdr(:,2))+1/n,max(bdr(:,2))-1/n,round(sqrt(n)));
[X,Y] = meshgrid(x,y);
pts = [X(:), Y(:)];
pts = pts + 1/n*randn(size(pts));
[keep,rem] = inpolygon(X(:),Y(:), bdr(:,1), bdr(:,2));
pts = pts(keep&~rem,:);
plot(pts(:,1), pts(:,2),'.','markersize',15)


%% Create grid
G = clippedPebi2D(pts, bdr);
Gc = createCVD(pts, bdr);

%% Plot
plotGrid(G)
axis equal off tight
figure()
plotGrid(Gc)
axis equal off tight
%% Save grids
save('showInnerOuterCellsSheep','G','Gc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Elephant
%% Load domain

mat = load('elephant.mat');
bdr = mat.P;
bdr = bdr/norm(bdr,'inf');

plot(bdr(:,1), bdr(:,2));
hold on
%% Create seeds.
rng(3)
n = 300;
x = linspace(min(bdr(:,1))+1/n,max(bdr(:,1))-1/n,round(sqrt(n)));
y = linspace(min(bdr(:,2))+1/n,max(bdr(:,2))-1/n,round(sqrt(n)));
[X,Y] = meshgrid(x,y);
pts = [X(:), Y(:)];
pts = pts + 1/n*randn(size(pts));
ptsx = (max(x)-min(x))*rand(n,1) + min(x); 
ptsy = (max(y)-min(y))*rand(n,1) + min(y); 
pts = [ptsx, ptsy];
[keep,rem] = inpolygon(ptsx ,ptsy,bdr(:,1), bdr(:,2));
pts = pts(keep&~rem,:);
plot(pts(:,1), pts(:,2),'.','markersize',15)


%% Create grid
G = clippedPebi2D(pts, bdr);
[Gc,optPts,f,g] = createCVD(pts, bdr);

%% Plot
plotGrid(G)
axis equal off tight
figure()
plotGrid(Gc)
axis equal off tight
%% Save grids
save('showInnerOuterCellsElephant')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot inner and outer cells
bndrFace = any(Gc.faces.neighbors==0,2); 
outer = [Gc.faces.neighbors(Gc.faces.neighbors(:,1)==0,2); ...
        Gc.faces.neighbors(Gc.faces.neighbors(:,2)==0,1)];
outer = unique(outer);

col = get(gca,'ColorOrder');
figure()
plotGrid(Gc,'facecolor','none')
plotGrid(Gc,outer,'facecolor','y')

for i = 1:numel(outer)
  c = col(mod(i, size(col,1))+1,:);
  pltF = any(Gc.faces.neighbors==outer(i),2)&any(Gc.faces.neighbors==0,2);
  plotFaces(Gc,pltF,'linewidth',2,'edgecolor',c)
end
axis equal off tight

