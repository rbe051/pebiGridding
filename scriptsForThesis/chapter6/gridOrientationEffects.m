clear;close all
addpath ../../voronoi3D/
addpath ../../../distmesh/
addpath ../
addpath ../../../vem/mat/
set(0,'defaulttextinterpreter','none')
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

fGs = x/21;
dt = x/21;
rho = @(p) 0.8*fGs*ones(size(p,1),1);

%% Create fault

% Add wells an faults as fixed points
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;

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


fDt = {fDt2};
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
%G = restrictedVoronoiDiagram(pts,bdrDT);
G = voronoi3D(pts,bdr);
G = computeGeometry(G);

%% CartGrid

rx = 0:dt:x;
ry = 0:dt:y;
rz = 0:dt:z;
nx = numel(rx);
ny = numel(ry);
nz = numel(rz);

slope = faultHeight2([0,1]) - faultHeight2([0,0]);
bound = [0,0; 0,y; x,y;x,0];

Gl = cartGrid([nx,ny,nz], [x,y,z]);
Gc = Gl;
Gc.nodes.coords(:,1) = Gl.nodes.coords(:,1)+slope*Gl.nodes.coords(:,1).*(x-Gl.nodes.coords(:,1)).*(Gl.nodes.coords(:,3) - z/2)/(x^2/4);
Gc = computeGeometry(Gc);
plotGrid(Gc)
axis equal

%% set bc
tol = 1e-6;
val = 1e8;
leftBc = Gc.faces.centroids(:,1)<tol;
rightBc = Gc.faces.centroids(:,1)>x -tol;
bcc=addBC([],find(leftBc),'pressure',val);
bcc=addBC(bcc,find(rightBc),'pressure',0);


leftB = G.faces.centroids(:,1)<tol;
rightB = G.faces.centroids(:,1)>x -tol;
bc=addBC([],find(leftB),'pressure',val);
bc=addBC(bc,find(rightB),'pressure',0);

%% Init fluid and rock
rockc.perm = milli * darcy * ones(Gc.cells.num,3);         % corner point
rockc.poro =  0.001 * ones(Gc.cells.num,1);                

rock.perm = milli * darcy * ones(G.cells.num,3);         %  unstructured
rock.poro =  0.001 * ones(G.cells.num,1);


fluid = initSingleFluid('mu' , 1*centi*poise     , ...
                        'rho', 1014*kilogram/meter^3);

Tc   = computeTrans(Gc, rockc);
T    = computeTrans(G, rock);

stateInitc = initState(Gc,{},0);
stateInit = initState(G,{},0);  

statec = incompTPFA(stateInitc,Gc,Tc,fluid,'bc',bcc);
state = incompTPFA(stateInit,G,T,fluid,'bc',bc);

%% plot
close all
figure(1); hold on
plotCellData(Gc,statec.pressure)
axis equal tight off
view(0,0)
figure(2); hold on
plotCellData(G,state.pressure)
axis equal tight off
view(0,0)

l = faultHeight2([0,0;0,z]);
l = [l, [-0.01;-0.01],[0;z]];
plot3(l(:,1), l(:,2), l(:,3),'r')
figure(1)
plot3(l(:,1), l(:,2), l(:,3),'r')

figure()
colorbar()
set(gca, 'CLim', [0, val]);
%% save
figure(1)
pause(0.1)
print('../../../../master/thesis/fig/ch06/gridOrientationEffectsCornerPtn','-depsc')
figure(2)
pause(0.1)
print('../../../../master/thesis/fig/ch06/gridOrientationEffectsUnstrctr','-depsc')
figure(3)
pause(0.1)
print('../../../../master/thesis/fig/ch06/inkscape/colorBarForGridOrientationEff','-dsvg')
