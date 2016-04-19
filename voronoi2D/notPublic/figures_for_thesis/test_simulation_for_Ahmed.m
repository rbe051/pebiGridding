close all; clear;
%% Grid parameters
xmax = 20;                              % Set grid dimentions
ymax = 20; 
gridSize = xmax*1/10;                   % distmesh

faultLine = {[16, 5; 3,10.1]};
wellLine = {[5,5], [15.0,15.0]};        % Set source center

mlqtMax = 2;                            % Set number of reminement levels
wellGridSize = 0.75/2^mlqtMax;          % Set relative well grid size
mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                        % Size around wells to be refined
                              
%% Set simulation parameters
T      = 120*second();    % End time
dT     = T/120;           % Time steps
dTplot = 1:1:T;           % Plot at these time steps


%%Generate grid

G = compositePebiGrid(gridSize, [xmax, ymax], 'wellLines', wellLine, ...
                     'wellGridFactor', wellGridSize, ...
                     'mlqtMaxLevel', 2, 'mlqtLevelSteps', mlqtSizes,...
                     'faultLines', faultLine);
                 
%Set internal boundary
G = makeInternalBoundary(G, find(G.faces.tag));
% Create 3d grid
G3D = makeLayeredGrid(G, 5);
G3D = computeGeometry(G3D);

%% Set fluid and rock properties
gravity reset off 

fluid = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
                        'rho', [1000, 700]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

rock.poro = ones(G3D.cells.num,1)*0.15;
rock.perm = ones([G3D.cells.num,1])*100*milli*darcy;


%% Create wells
wells = find(G3D.cells.tag);
pv = sum(poreVolume(G3D, rock));
W = verticalWell([], G3D, rock, wells(1),[], 'radius', wellGridSize/10,...
                'type', 'rate', 'val', 0.01*pv, 'Comp_i', [1,0], 'name','$I$');
W = verticalWell(W, G3D, rock, wells(2),[], 'radius', wellGridSize/10,...
                'type','bhp', 'val', 300*darcy(),'Comp_i', [0,1], 'name','$P$' );  
            
            

%% Solve
state = initState(G3D, W, 0, [0.0,1]);
trans = computeTrans(G3D,rock);
state = incompTPFA(state, G3D, trans, fluid, 'wells', W);
% Prepare plotting of saturations
clf;
hold on
plotGrid(G3D, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
plotWell(G3D, W, 'height', 10, 'color', 'c');
axis off equal, view([-120,30]), colormap(flipud(jet))

colorbar; hs = []; ha=[]; zoom(1.3);

% Start the main loop
t  = 0;  plotNo = 1;
while t < T,
   state = implicitTransport(state, G3D, dT, rock, fluid, 'wells', W);

   % Check for inconsistent saturations
   assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);
   
   % Update solution of pressure equation.
   %state = incompMimetic(state, G3D, S, fluid, 'wells', W);
   state = incompTPFA(state, G3D, trans, fluid, 'wells', W);
   
   % Increase time and continue if we do not want to plot saturations
   t = t + dT;

   
   if ( t + dT <= dTplot(plotNo)), continue, end
   delete([hs, ha])
   hs = plotCellData(G3D, state.s(:,1), find(state.s(:,1) >= 0.0));
   ha = annotation('textbox', [0.1714 0.8214 0.5000 0.1000], 'LineStyle', 'none', ...
                   'String', ['Water saturation at ', ...
                              num2str(convertTo(t,second)), ' s']);

   fig = gcf();
   set(findall(fig,'-property','FontSize'),'FontSize',14) 
   view(0, 90), drawnow, caxis([0 1])
  
   plotNo = plotNo+1;
end