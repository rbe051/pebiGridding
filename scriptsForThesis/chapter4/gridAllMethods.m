clear; close all
%addpath ../../voronoi2D/
%addpath ../
%addpath ../../voronoi3D/
%addpath ../../../vem/mat/
%addpath ../../../vem/mat/VEM2D
%addpath ../../../distmesh/

%% Set gridding parameters
fGs = 0.05;
wGs = 0.04;
rGs = 0.05;
%% Set fault and well
f ={[0.2,0.2; 0.5,0.6; 0.8,0.8],...
    [0.5,0.6; 0.8,0.5]};
  
w = {[0.2,0.9; 0.5,0.2], ...
     [0.2,0.7; 0.7,0.8]};
plotLinePath(f)
plotLinePath(w)

[fs,fCut,fwCut] = splitAtInt(f,w);
[ws,wCut,wfCut] = splitAtInt(w,f);


%% create well and fault sites
bisectPnt = (fGs.^2 - (0.6*fGs).^2 ...
            + (0.6*fGs).^2) ./(2*fGs);
faultOffset = sqrt((0.6*fGs).^2 - bisectPnt.^2);
sePtn = [wfCut==2|wfCut==3, wfCut==1|wfCut==3];
sePtn = (1.0+faultOffset/wGs)*sePtn;

[wSites,wGs] = createWellGridPoints(ws, wGs,'sePtn',sePtn);

F = createFaultGridPoints(fs,fGs,'fwCut',fwCut,'fCut',fCut);

%% Create cartesian Grid
bnd = [0,0;0,1;1,1;1,0];
[X,Y]  = meshgrid(0:fGs*0.8:1);
rSites = [X(:),Y(:)];
rSites = removeConflictPoints2(rSites, [F.f.pts;wSites],[F.f.Gs;wGs]);
pts = [F.f.pts;wSites;rSites];
Gc = clippedPebi2D(pts, bnd);


%% Create distmesh Grid
x = [1,1];
rectangle = [0,0; x(1),x(1)];   
fd = @(p) drectangle(p, 0, x(1), 0, x(2));
% Set fixed points
corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];
fixedPts = [F.f.pts;wSites;corners];
% Set density function
hd = @(p) ones(size(p,1),1);

[Pts,t,sorting] = distmesh2d(fd, hd, rGs, rectangle, fixedPts);

Gd = clippedPebi2D(Pts,bnd);

%% Create CVD Grid
n = 560;%round((norm(x)/rGs)^2);
rSi = rand(n,2);
wAndF = [F.f.pts;wSites];
nf = size(wAndF,1);
[G,optPts,func,g] = createCVD(rSi,bnd,'fixedPts',wAndF);
newOpt = removeConflictPoints2(optPts(nf+1:end,:),F.f.pts,F.f.Gs);
newOpt = [wAndF;newOpt];
Gcvd = clippedPebi2D(newOpt,bnd);

%% Plot Grid
close all
color = get(gca,'colororder');
figure(1)
plotGrid(Gc,'facecolor','none')
plotLinePath(w,'--','color',color(1,:),'linewidth',1)
plotLinePath(f,'--','color',color(2,:),'linewidth',1)
axis off equal tight
figure(2)
plotGrid(Gd,'facecolor','none')
plotLinePath(w,'--','color',color(1,:),'linewidth',1);
plotLinePath(f,'--','color',color(2,:),'linewidth',1)
axis off equal tight
figure(3)
plotGrid(Gcvd,'facecolor','none')
plotLinePath(w,'--','color',color(1,:),'linewidth',1);
plotLinePath(f,'--','color',color(2,:),'linewidth',1)
axis off equal tight

%% Save grids
figure(1)
print('../../../../master/thesis/fig/ch04/compareResPtsMethodCart','-depsc')
pause(0.1)
figure(2)
print('../../../../master/thesis/fig/ch04/compareResPtsMethodDist','-depsc')
pause(0.1)
figure(3)
print('../../../../master/thesis/fig/ch04/compareResPtsMethodCVD','-depsc')
