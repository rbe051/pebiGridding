clear;close all
addpath ../../voronoi3D/
addpath ../../../distmesh/

%% set boundary
x = 10;
y = 10;
z = 5;
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
rho = @(p) 0.8*fGs*ones(size(p,1),1);

%% Create fault
rectangle = [min(bdr(:,1:2)); max(bdr(:,1:2))-0];
% Add wells an faults as fixed points
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
fixedPts = [rectangle; rectangle(swap1)';rectangle(swap2)'];
hd = @(p) ones(size(p,1),1);
fd = @(p) drectangle(p, rectangle(1),rectangle(2), rectangle(3),rectangle(4));

%fault 1
[Pts1,t1] = distmesh2d(fd, hd, fGs, rectangle, fixedPts);
fDt1.ConnectivityList = t1;

faultHeight1 = @(p) y/4 + y/2*sin(pi*p(:,1)/(2*x));
fDt1.Points = zeros(size(Pts1,1),3);
fDt1.Points = [Pts1(:,1), faultHeight1(Pts1),Pts1(:,2)];

% fault 2
[Pts2,t2] = distmesh2d(fd, hd, fGs, rectangle, fixedPts);
fDt2.Points = Pts2;
fDt2.ConnectivityList = t2;
faultHeight2 = @(p) 3*z/4 - z/2*sin(pi*p(:,1)/(2*x));
fDt2.Points = [fDt2.Points, faultHeight2(fDt2.Points)];


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
close all
figure()
plotGrid(Gd)
axis equal
figure()
c = Gd.cells.centroids(:,3)>faultHeight2(Gd.cells.centroids(:,1:2));
plotGrid(Gd,c)
axis equal