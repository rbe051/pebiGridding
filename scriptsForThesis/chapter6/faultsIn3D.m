clear;close all
%% set boundary
x = 10;
y = 3;
z = 2;
bdr   = [ 0, 0, 0;  ...
          x, 0, 0;  ...
          x, y, 0;  ...
          0, y, 0;  ...
          0, 0, z;  ...
          x, 0, z;  ...
          x, y, z;  ...
          0, y, z];

%% Set gridding parameters
n = 1000;
fGs = x/40;
rho = @(p) fGs*0.8*ones(size(p,1),1);

%% Create fault
rectangle1 = [min(bdr(:,1:2)); max(bdr(:,1:2))-0];
% Add wells an faults as fixed points
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
fixedPts = [rectangle1; rectangle1(swap1)';rectangle1(swap2)'];
hd = @(p) ones(size(p,1),1);
fd = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3),rectangle1(4));

[Pts1,t1] = distmesh2d(fd, hd, fGs, rectangle1, fixedPts);
fDt1.Points = Pts1;
fDt1.ConnectivityList = t1;
faultHeight1 = @(p) z/4 + z/2*sin(pi*p(:,1)/(2*x));
fDt1.Points = [fDt1.Points, faultHeight1(fDt1.Points)];

rectangle2 =[min(bdr(:,1)),4;max(bdr(:,1)),6];
fixedPts2 = [rectangle2; rectangle2(swap1)';rectangle2(swap2)'];
fd2 = @(p) drectangle(p, rectangle2(1),rectangle2(2), rectangle2(3),rectangle2(4));
[Pts2,t2] = distmesh2d(fd2, hd, fGs, rectangle2, fixedPts2);
fDt2.Points = Pts2;
fDt2.ConnectivityList = t2;
faultHeight2 = @(p) (rectangle2(2,2) - p(:,2))*max(bdr(:,3))/(rectangle2(2,2)-rectangle2(1,2));
fDt2.Points = [fDt2.Points, faultHeight2(fDt2.Points)];


fDt = {fDt1,fDt2};
F = createFaultGridPoints3D(fDt,rho);

%% Create reservoir sites

%rSites = rand(n,3);
xmax = max(bdr(:,1))+0.01; xmin = min(bdr(:,1))-0.01;
ymax = max(bdr(:,2))+0.01; ymin = min(bdr(:,2))-0.01;
zmax = max(bdr(:,3))+0.01; zmin = min(bdr(:,3))-0.01;

dt = 0.5*(1/n^(1/3));
x = xmin:0.5*dt:xmax;
y = ymin:dt:ymax;
z = zmin:dt:zmax;
[X,Y,Z] = ndgrid(x,y,z);
rSites = [X(:), Y(:), Z(:)];
rSites = bsxfun(@times, rSites,max(bdr) - min(bdr));
rSites = bsxfun(@minus, rSites, min(bdr));

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
%Gd = restrictedVoronoiDiagram(pts,bdrDT);
Gd = voronoi3D(pts,bdr);
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