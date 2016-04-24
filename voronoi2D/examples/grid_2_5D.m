%% Example
% This script contains an example of how one can create a 2.5D grid.
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.


%% Set path to voronoi2D
%addpath ../
%close all; clear

%% Set well and fault paths
% We start by creating wells and faults. We create three vertical wells and
% two faults. We will first grid one layer in the 2.5D grid and then
% extrude this to obtain the final grid. 

[nx,ny,nz] = deal(20,10,5);    % Number of elements in x,y,z directions
[Lx,Ly,Lz] = deal(400,200,50);  % Length of domain in x,y,z directions

% Create the wells and faults.
fault = {[40,60;360,180], [40,160;250,40]};
well  = {[40,100], [350, 150], [160,140]};

%% Plot faults and well paths
% and plot them
figure(); hold on
plotLinePath(well,'.','color','blue','markersize',15);
plotLinePath(fault,'color','red');
axis equal tight
axis ([0,Lx,0,Ly])
title('well & fault paths')
drawnow


%% Setting gridding parameters:
% Before we call the gridding fuctions we set some parameters.
gS  = norm([Lx,Ly])/norm([nx,ny]); % set grid size
wGf = 0.5;  % The relative size of well cells compared to gS
fGf = 0.5;  % The relative size of fault cells compared to gS
nRs = 1;    % Number of refinement steps towards the well

%% Create Grid
% We can now create the composite Pebi grid:
Gc = compositePebiGrid(gS, [Lx, Ly], 'wellLines', well, ...
                               'wellGridFactor', wGf,...
                               'mlqtMaxLevel', nRs, 'faultLines', fault);
%% Plot composite grid
%And plot it
figure(); hold on
plotGrid(Gc, 'facecolor','none')
axis equal tight
title('compositePebiGrid(...)')
drawnow

%% Set pebiGrid Parameters:
% We now use the other wrapper function to create a PEBI-grid using
% distmesh:
eps = sqrt(gS)*8; % This parameter defines the refinement around 
                  % the wells. The cell size increase exponentialy from the
                  % wells as exp(distance from well/eps)

%% Generate grid
% distmesh will most likely not converge in the maximum number of iterations.
% This is usually not a problem, since the grid most likely is good before
% the convergence requirement is met. 

Gdist = pebiGrid(gS, [Lx, Ly], 'wellLines', well, ...
                'wellGridFactor',wGf/2, 'wellRefinement', true, ...
                'epsilon',eps, 'faultlines', fault);
%% Plot pebiGrid
figure(); hold on
plotGrid(Gdist,'facecolor','none')
axis equal tight
title('pebiGrid(...)')

%% Create 2.5D grid
% We create a 3dGrid by extrapolate the 2D grid. 
Gc3D    = makeLayeredGrid(Gc, nz);
Gc3D.nodes.coords(:,3) = Gc3D.nodes.coords(:,3)*Lz/nz;
Gc3D = computeGeometry(Gc3D);

Gdist3D = makeLayeredGrid(Gdist, nz);
Gdist3D.nodes.coords(:,3) = Gdist3D.nodes.coords(:,3)*Lz/nz;
Gdist3D = computeGeometry(Gdist3D);

%% plot 3D grids
figure(); hold on
plotGrid(Gc3D);
axis equal tight
title('2.5D grid')
view(30,60)
        
figure(); hold on
plotGrid(Gdist3D);
axis equal tight
title('2.5D grid')
view(30,60)        