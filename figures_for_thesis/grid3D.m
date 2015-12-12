%% Grid with single source centered at 0.5,0.5
close all; clear;
addpath('../', '../mrstTweak/')

typeOfGrid = 'coarseCart';

%% grid parameters
xmax = 10;                              % Set grid dimentions
ymax = 10; 
gridSize = 10*1/15;                       % Set size of grid cells

faultLine = {};
wellLine = {[2,2], [8.0,8.0]};                % Set source center

mlqtMax = 2;                        % Set number of reminement levels
wellGridSize = 0.5/2^mlqtMax;       % Set relative well grid size
mlqtSizes = 2*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                    % Set distance around wells to be
                                    % refined.

wellEps = sqrt(10)*1/4;                         % Size around wells to be refined
                                       % (For unstructured grid)
%% Set simulation parameters
T      = 60*second();
dT     = T/100;
dTplot = T/4;


%% Generate grid
if strcmp(typeOfGrid, 'composite')
    % Create semi-structured grid
    G = compositeGridPEBI(gridSize, [xmax, ymax], 'wellLines', wellLine, ...
                         'wellGridFactor', wellGridSize, ...
                         'mlqtMaxLevel', 2, 'mlqtLevelSteps', mlqtSizes,...
                         'padding', 1);
elseif strcmp(typeOfGrid, 'distmesh')
    %Create fully unstructured grid
    G = compositeGridPEBIdistmesh(gridSize, [xmax, ymax], 'wellLines', wellLine, ...
                                 'wellGridFactor', wellGridSize, ...
                                 'wellRefDist' ,wellEps, 'padding',0);
elseif strcmp(typeOfGrid,'coarseCart')
    nx = ceil(xmax/gridSize);
    ny = ceil(ymax/gridSize);
    G = cartGrid([nx,ny],[xmax,ymax]);
    
end

% Create 3d grid
G3D    = makeLayeredGrid(G   ,5);
G3D = computeGeometry(G3D);


%% Set fluid properties
gravity reset off 

fluid = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
                        'rho', [1000, 700]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%% Set rock
rock.poro = ones(G3D.cells.num,1)*0.2;
rock.perm = ones([G3D.cells.num,1])*100*milli*darcy;


%% Create Wells
wells = find(G3D.cells.isWell);
pv = sum(poreVolume(G3D, rock));
W = verticalWell([], G3D, rock, wells(1),[], 'radius', wellGridSize/10,...
                'type', 'rate', 'val', 0.1*pv, 'Comp_i', [1,0], 'name','$I$');
W = verticalWell(W, G3D, rock, wells(2),[], 'radius', wellGridSize/10,...
                'type','bhp', 'val', 300*darcy(), 'Comp_i', [0,1], 'name','P' ) ;

%% Reservoir states
state = initState(G3D, W, 0, [0.1,1]);

%% Solve
S = computeMimeticIP(G3D, rock);
state = incompMimetic(state, G3D, S, fluid, 'wells', W);

% Prepare plotting of saturations
clf
plotGrid(G3D, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
plotWell(G3D, W, 'height', 2, 'color', 'c');
axis off equal, view([-120,30]), colormap(flipud(jet))
colorbar; hs = []; ha=[]; zoom(1.3);
 
% Start the main loop
t  = 0;  plotNo = 0;
while t < T,
   state = implicitTransport(state, G3D, dT, rock, fluid, 'wells', W);

   % Check for inconsistent saturations
   assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);

   % Update solution of pressure equation.
   state = solveIncompFlow(state, G3D, S, fluid, 'wells', W);

    % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t), continue, end

   % Plot saturation
   delete([hs, ha])
   hs = plotCellData(G3D, state.s(:,1), find(state.s(:,1) > 0.01));
   ha = annotation('textbox', [0.1714 0.8214 0.5000 0.1000], 'LineStyle', 'none', ...
                   'String', ['Water saturation at ', ...
                              num2str(convertTo(t,second)), ' seconds']);
   view(-105, 40), drawnow, caxis([0 1])
   plotNo = plotNo+1;
   pause()
end


