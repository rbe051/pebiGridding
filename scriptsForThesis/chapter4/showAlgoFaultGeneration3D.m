clear; close all
%addpath ../../voronoi3D/

%% set boundary
x = 1;
y = 1;
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
fGs = x/40;
dt = x/20;
rho = @(p) fGs*(1+4*p(:,1));

%% Create Fault
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
rectangle1 = [min(bdr(:, [1,3])); max(bdr(:,[1,3]))];
fixedPts1 = [rectangle1; rectangle1(swap1)';rectangle1(swap2)'];
hd = @(p) rho(p)/fGs;
fd1 = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3),rectangle1(4));

[Pts1,t1] = distmesh2d(fd1, hd, fGs, rectangle1, fixedPts1);

fDt1.ConnectivityList = t1;

faultHeight1 = @(p) y/3*ones(size(p,1),1) + 0.3*p(:,1)+0.2*p(:,2);
fDt1.Points = [Pts1(:,1), faultHeight1(Pts1), Pts1(:,2)];


fDt = {fDt1};
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
%pts = [X(:),Y(:),Z(:)];
bdrDT = delaunayTriangulation(bdr);
Gd = restrictedVoronoiDiagram(pts,bdrDT);
%Gd = voronoi3D(pts,bdr);
Gd = computeGeometry(Gd);

%% plot
color = get(gca,'ColorOrder');
close all

% Side view
figure(1); hold on
patch('Vertices', fDt1.Points, 'faces', fDt1.ConnectivityList,'facecolor','y')
axis equal;
% [X,Y,Z] = sphere(30);
% 
% for i = 1:size(F.c.CC,1)
%   p = F.c.CC(i,:);
%   R = F.c.R(i);
%   Xk = R*X + p(1);
%   Yk = R*Y + p(2);
%   Zk = R*Z + p(3);
%   surf(Xk,Yk,Zk,'edgecolor','none','facealpha',0.9);
% end
plot3(F.f.pts(:,1), F.f.pts(:,2),F.f.pts(:,3),'.','color',color(2,:),'markersize',7)
view(0,0)
axis([0,1,0,1,0,1])
axis equal tight off

% Top view
figure(2); hold on
patch('Vertices', fDt1.Points, 'faces', fDt1.ConnectivityList,'linewidth',2)
plot3(F.f.pts(:,1), F.f.pts(:,2),F.f.pts(:,3),'.','color',color(2,:),'markersize',7)
view(19,-79)
axis([0,1,0,1,0,1])
axis equal tight off


% Grid view
figure(3); hold on
c = Gd.cells.centroids(:,2)>faultHeight1(Gd.cells.centroids(:,[1,3]));
plotGrid(Gd,c)
 % Plot box
 x = [0,1,1,0,0];
y = [0,0,1,1,0];
z = [0,0,0,0,0];

  plot3(x,y,z,'k')
  plot3(x,y,z+1,'k')
  for j = 1:numel(x)
    plot3([x(j),x(j)],[y(j),y(j)],[0,1],'k')
  end


view(60,30);
light('Position',[-1 -1 -1],'Style','infinite')
axis equal tight off

%% Save plots
figure(1)
print('../../../../master/thesis/fig/ch04/showAlgo3DFaultSide','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch04/showAlgo3DFaultTop','-depsc')
figure(3)
print('../../../../master/thesis/fig/ch04/showAlgo3DFaultGrid','-depsc')

