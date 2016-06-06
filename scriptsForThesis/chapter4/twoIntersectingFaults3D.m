clear;close all
%addpath ../../voronoi3D/
%addpath ../../../distmesh/

%% set boundary
x = 5;
y = 5;
z = 1;
bdr   = [ 0, 0, 0;  ...
          x, 0, 0;  ...
          x, y, 0;  ...
          0, y, 0;  ...
          0, 0, z;  ...
          x, 0, z;  ...
          x, y, z;  ...
          0, y, z];

%% Set gridding parameters

fGs = x/20;
dt = x/20;
rho = @(p) fGs*ones(size(p,1),1);

%% Create fault
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
rectangle1 = [min(bdr(:, [1,3])); max(bdr(:,[1,3]))];
fixedPts1 = [rectangle1; rectangle1(swap1)';rectangle1(swap2)'];
hd = @(p) 0.6*rho(p)/fGs;
fd1 = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3),rectangle1(4));

[Pts1,t1] = distmesh2d(fd1, hd, fGs, rectangle1, fixedPts1);
[Pts2,t2] = distmesh2d(fd1, hd, fGs, rectangle1, fixedPts1);
fDt1.ConnectivityList = t1;
fDt2.ConnectivityList = t2;

faultHeight1 = @(p) y/3*ones(size(p,1),1) + 0.3*p(:,1);
faultHeight2 = @(p) 2*y/3*ones(size(p,1),1) - x/2*sin(pi*p(:,1)/(2*x));

fDt1.Points = [Pts1(:,1), faultHeight1(Pts1), Pts1(:,2)];
fDt2.Points = [Pts2(:,1), faultHeight2(Pts2), Pts2(:,2)];


fDt = {fDt1,fDt2};
F = createFaultGridPoints3D(fDt,rho);


%% Create reservoir sites

%rSites = rand(n,3);
xmax = max(bdr(:,1))-0.001; xmin = min(bdr(:,1))+0.001;
ymax = max(bdr(:,2))-0.001; ymin = min(bdr(:,2))+0.001;
zmax = max(bdr(:,3))-0.001; zmin = min(bdr(:,3))+0.001;


x = xmin:dt:xmax;
y = ymin:dt:ymax;
z = zmin:dt:zmax;

[X,Y,Z] = ndgrid(x,y,z);
rSites = [X(:), Y(:), Z(:)];

nR = size(rSites,1);
rGs = zeros(nR,1);
rPri = zeros(nR,1);

[rSites,removed] = faultSufCond(rSites,F.c.CC,F.c.R);
rGs = rGs(~removed);
rPri = rPri(~removed);

for i = 1:size(F.l.fPos)-2
  f1 = F.l.f(F.l.fPos(i):F.l.fPos(i+1)-1);
  c2 = F.l.c(F.l.cPos(i+1):F.l.cPos(i+2)-1);
  [~,removed] = faultSufCond(F.f.pts(f1,:),F.c.CC(c2,:),F.c.R(c2));
  F.f.pts(f1(removed),:) = []; 
  LIA = ismember(F.l.f,find(removed));
  sub = cumsum(LIA);
  F.l.f = F.l.f - sub;
  F.l.f(LIA) = [];
  F.l.fPos(2:end) = F.l.fPos(2:end) - sub(F.l.fPos(2:end)-1);
end
pts = [F.f.pts;rSites];
gs  = [F.f.Gs;rGs];
pri = [F.f.pri;rPri];
%[pts,removed] = removeConflictPoints3(pts,4*gs, pri );

bdrDT = delaunayTriangulation(bdr);
Gd = restrictedVoronoiDiagram(pts,bdrDT);
%Gd = voronoi3D(pts,bdr);
Gd = computeGeometry(Gd);
%% plot grid
color = get(gca,'ColorOrder');
close all
% top view triangulation
figure(1); hold on
for i = 1:numel(fDt)
  patch('Vertices', fDt{i}.Points, 'faces', fDt{i}.ConnectivityList,'linewidth',2)
  axis equal;
  f = F.l.f(F.l.fPos(i):F.l.fPos(i+1)-1);
  plot3(F.f.pts(f,1), F.f.pts(f,2),F.f.pts(f,3),'.','color',color(1+i,:),'markersize',7)
  view(0,-90)
  axis([0,1,0,1,0,1])
  axis equal tight off
end

%%
figure(2)
plotGrid(Gd)

view(38,45);
light('Position',[6, -6 -5],'Style','local')
light('Position',[6, -6 -5],'Style','local')
axis equal tight off

%% Save plots
figure(1)
print('../../../../master/thesis/fig/ch04/inter2Faults3dTop','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch04/inter2Faults3dGrid','-depsc')


