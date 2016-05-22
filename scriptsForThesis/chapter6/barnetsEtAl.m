clear; close all

%% sett geometry
gridBranetsEtAl;

%% calculate fault intersections
[flt, fCut, fwCut] = splitAtInt(flt, {});

%% Create fault seeds
fGs = max(max(bdr))/70;
F = createFaultGridPoints(flt, fGs,'fCut',fCut,'fwCut', fwCut);
figure(); hold on
plot(F.f.pts(:,1), F.f.pts(:,2),'.','markersize',15);

%% Create reservoir seeds
% n = 2000;
% x = linspace(min(bdr(:,1))+1/n,max(bdr(:,1))-1/n,round(sqrt(n)));
% y = linspace(min(bdr(:,2))+1/n,max(bdr(:,2))-1/n,round(sqrt(n)));
% [X,Y] = meshgrid(x,y);
% rPts = [X(:), Y(:)];
% rPts = rPts + 1/n*randn(size(rPts));
% [keep,rem] = inpolygon(X(:),Y(:), bdr(:,1), bdr(:,2));
% rPts = rPts(keep&~rem,:);
% 
% rPts = removeConflictPoints2(rPts,F.f.pts, F.f.Gs);
% plot(rPts(:,1), rPts(:,2),'.','markersize',15)
% Create Reservoir grid points
% set domain function

rectangle = [min(bdr); max(bdr)];   

% Add wells an faults as fixed points
fixedPts = [F.f.pts;bdr];
uni = @(p,varargin) 2*ones(size(p,1),1);
[Pts,~,sorting] = distmesh2d(@signDistFuncBarn ,uni, fGs, rectangle, fixedPts,'bdr',bdr);

% Distmesh change the order of all points. We undo this sorting.
isFault = false(size(Pts,1),1); isFault(1:size(F.f.pts,1)) = true;
isFault = isFault(sorting);
[~,If]   = sort(sorting(isFault));
isRes   = ~isFault;

fPts = Pts(isFault,:);
fPts = fPts(If,:);
rPts = Pts(isRes,:);
  

%% Create grid
pts = [F.f.pts;rPts];
G = clippedPebi2D(pts,bdr);

%% Plot grid
plotGrid(G)
axis equal off tight
%plot(F.f.pts(:,1), F.f.pts(:,2),'.','markersize',15)
%plotLinePath(flt)