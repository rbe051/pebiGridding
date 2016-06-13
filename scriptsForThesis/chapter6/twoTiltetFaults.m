clear;close all
addpath ../../voronoi3D/
addpath ../../../distmesh/

%% set boundary
x = 20;
y = 10;
z = 10;
bdr   = [ 0, 0, 0;  ...
          x, 0, 0;  ...
          x, y, 0;  ...
          0, y, 0;  ...
          0, 0, z;  ...
          x, 0, z;  ...
          x, y, z;  ...
          0, y, z];

%% Set gridding parameters

fGs = x/40;
dt = x/40;
rho = @(p) 0.8*fGs*ones(size(p,1),1);

%% Create fault

% Add wells an faults as fixed points
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
% fault 1
rectangle1 = [min(bdr(:,[1,2])); 5/9*x,max(bdr(:,2))];
fixedPts1 = [rectangle1; rectangle1(swap1)';rectangle1(swap2)'];
hd = @(p) ones(size(p,1),1);
fd1 = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3),rectangle1(4));

[Pts1,t1] = distmesh2d(fd1, hd, fGs, rectangle1, fixedPts1);
fDt1.ConnectivityList = t1;

faultHeight1 = @(p) 2*z/3*ones(size(p,1),1);
fDt1.Points = zeros(size(Pts1,1),3);
fDt1.Points = [Pts1(:,1), Pts1(:,2),faultHeight1(Pts1)];


% fault 2
rectangle2 = [min(bdr(:, [2,3])); max(bdr(:,[2,3]))];
fixedPts2 = [rectangle2; rectangle2(swap1)';rectangle2(swap2)'];
hd = @(p) ones(size(p,1),1);
fd2 = @(p) drectangle(p, rectangle2(1),rectangle2(2), rectangle2(3),rectangle2(4));

[Pts2,t2] = distmesh2d(fd2, hd, fGs, rectangle2, fixedPts2);
fDt2.ConnectivityList = t2;

faultHeight2 = @(p) x/3 + 1/3*x*p(:,2)/z;
fDt2.Points = zeros(size(Pts2,1),3);
fDt2.Points = [faultHeight2(Pts2),Pts2(:,1), Pts2(:,2)];

% fault 3
rectangle3 = [4/9*x,min(bdr(:,2)); max(bdr(:,[1,2]))];
fixedPts3 = [rectangle3; rectangle3(swap1)';rectangle3(swap2)'];
hd = @(p) ones(size(p,1),1);
fd3 = @(p) drectangle(p, rectangle3(1),rectangle3(2), rectangle3(3),rectangle3(4));

[Pts3,t3] = distmesh2d(fd3, hd, fGs, rectangle3, fixedPts3);
fDt3.ConnectivityList = t3;

faultHeight3 = @(p) z/3*ones(size(p,1),1);
fDt3.Points = zeros(size(Pts3,1),3);
fDt3.Points = [Pts3(:,1), Pts3(:,2),faultHeight3(Pts3)];


fDt = {fDt2};%,fDt3,fDt1};
F = createFaultGridPoints3D(fDt,rho);

%% Create reservoir sites

%rSites = rand(n,3);
xmax = max(bdr(:,1))-0.001; xmin = min(bdr(:,1))+0.001;
ymax = max(bdr(:,2))-0.001; ymin = min(bdr(:,2))+0.001;
zmax = max(bdr(:,3))-0.001; zmin = min(bdr(:,3))+0.001;


rx = xmin:dt:xmax;
ry = ymin:dt:ymax;
rz = zmin:dt:zmax;
[X,Y,Z] = ndgrid(rx,ry,rz);
rSites = [X(:), Y(:), Z(:)];

nR = size(rSites,1);
rGs = zeros(nR,1);
rPri = zeros(nR,1);

[rSites,removedRes] = faultSufCond(rSites,F.c.CC,F.c.R);
rGs = rGs(~removedRes);
rPri = rPri(~removedRes);

figure(); hold on
for i = 1:size(F.l.fPos)-2
  f2 = F.l.f(F.l.fPos(i+1):F.l.fPos(end)-1);
  c1 = F.l.c(F.l.cPos(i):F.l.cPos(i+1)-1);
  [~,removed] = faultSufCond(F.f.pts(f2,:),F.c.CC(c1,:),F.c.R(c1));
  F.f.pts(f2(removed),:) = []; 
  LIA = ismember(F.l.f,find(removed)+F.l.fPos(i+1)-1);
  sub = cumsum(LIA);
  F.l.f = F.l.f - sub;
  F.l.f(LIA) = [];
  F.l.fPos(2:end) = F.l.fPos(2:end) - sub(F.l.fPos(2:end)-1);
end
pts = [F.f.pts;rSites];
gs  = [F.f.Gs;rGs];
pri = [F.f.pri;rPri];
%[pts,removed] = removeConflictPoints3(pts,4*gs, pri );
%pts = [X(:),Y(:),Z(:)];

bdrDT = delaunayTriangulation(bdr);
figure(); hold on
G = clippedPebi3D(pts,bdrDT);
%G = voronoi3D(pts,bdr);
G = computeGeometry(G);

%% CartGrid
ptsc = [X(:),Y(:),Z(:)];
ptsc = [ptsc(removedRes,:);ptsc(~removedRes,:)];
Gc = voronoi3D(ptsc,bdr);
Gc = computeGeometry(Gc);
%% shift geometry
cr1 = G.cells.centroids(:,3)<z/3;
cr1 = cr1& G.cells.centroids(:,1)>faultHeight2(G.cells.centroids(:,2:3));
cr2 = G.cells.centroids(:,3)>2*z/3;
cr2 = cr2& G.cells.centroids(:,1)<faultHeight2(G.cells.centroids(:,2:3));
cr = cr1|cr2;
Gr = removeCells(G,cr);

cr1c = Gc.cells.centroids(:,3)<z/3;
cr1c = cr1c& Gc.cells.centroids(:,1)>faultHeight2(Gc.cells.centroids(:,2:3));
cr2c = Gc.cells.centroids(:,3)>2*z/3;
cr2c = cr2c& Gc.cells.centroids(:,1)<faultHeight2(Gc.cells.centroids(:,2:3));
crc = cr1c|cr2c;
Grc = removeCells(Gc,crc);


%% plot grid
close all
figure()
shift = G.cells.centroids(:,1)>faultHeight2(G.cells.centroids(:,2:3));
shiftc = Gc.cells.centroids(:,1)>faultHeight2(Gc.cells.centroids(:,2:3));
c1 = G.cells.centroids(:,3) - 10/3*heaviside(shift-0.5)<2*z/9;
c2 = ~c1& G.cells.centroids(:,3) - 10/3*heaviside(shift-0.5)<4*z/9;
c3 = ~c1& ~c2;
c1c = Gc.cells.centroids(:,3) - 10/3*heaviside(shiftc-0.5)<2*z/9;
c2c = ~c1c& Gc.cells.centroids(:,3) - 10/3*heaviside(shiftc-0.5)<4*z/9;
c3c = ~c1c& ~c2c;

n1 = ceil(sum(c1).^(1/3));
p1 = gaussianField([n1,n1,n1],[0.05,0.2],3,10);
p1 = p1(:);
p1 = p1(1:sum(c1));
p1c = p1(1:sum(c1c));

n2 = ceil(sum(c1).^(1/3));
p2 = gaussianField([n2,n2,n2],[0.4,0.6],3,10);
p2 = p2(:);
p2 = p2(1:sum(c2));
p2c = p2(1:sum(c2c));

n3 = ceil(sum(c1).^(1/3));
p3 = gaussianField([n3,n3,n3],[0.2,0.4],3,10);
p3 = p3(:);
p3 = p3(1:sum(c3));
p3c = p3(1:sum(c3c));

p = [p1;p2;p3];
pc = [p1c;p2c;p3c];

c = [find(c1);find(c2);find(c3)];
cc = [find(c1c);find(c2c);find(c3c)];
[~,I] = sort(c);
[~,Ic] = sort(cc);
p = p(I);
pc = pc(Ic);

pc = [pc(1:sum(removedRes));p(size(F.f.pts,1)+1:end)];

plotCellData(Gr,p(~cr));
light('position',[25,0,0],'style','local')
light('position',[16,-5,0],'style','local')
light('position',[25,-5,5],'style','local')
light('position',[5,5,-10],'style','local')
axis equal off tight
view(30,25)
figure()
plotCellData(Grc,pc(~crc))
light('position',[25,0,0],'style','local')
light('position',[16,-5,0],'style','local')
light('position',[25,-5,5],'style','local')
light('position',[5,5,-10],'style','local')
axis equal off tight
view(30,25)

% Plot auxillary cells
figure()
plotGrid(G,~cr);
plotGrid(G,cr,'facecolor','r','facealpha',1)
light('position',[25,0,0],'style','local')
light('position',[16,-5,0],'style','local')
light('position',[25,-5,5],'style','local')
light('position',[5,5,-10],'style','local')
axis equal off tight
view(30,25)

%% Save
figure(1)
pause(0.1)
%print('../../../../master/thesis/fig/ch06/singleFault3DExact','-depsc')

figure(2)
pause(0.1)
%print('../../../../master/thesis/fig/ch06/singleFault3DCart','-depsc')

figure(3)
pause(0.1)
%print('../../../../master/thesis/fig/ch06/singleFault3DAux','-depsc')


