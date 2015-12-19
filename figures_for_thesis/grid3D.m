%% Grid with single source centered at 0.5,0.5
close all; clear;
addpath('../', '../mrstTweak/')
mrstModule add mimetic
mrstModule add distmesh

typeOfGrid = 'distmesh';
fileFormat = 'pdfNoVector';

%% grid parameters
xmax = 20;                              % Set grid dimentions
ymax = 20; 
gridSize = xmax*1/7.5;                       % Set size of grid cells

faultLine = {[16, 5; 3,10.1]};
wellLine = {[5,5], [15.0,15.0]};                % Set source center

mlqtMax = 2;                        % Set number of reminement levels
wellGridSize = 0.75/2^mlqtMax;       % Set relative well grid size
mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                    % Set distance around wells to be
                                    % refined.

wellEps = sqrt(xmax)*1/2;                         % Size around wells to be refined
                                       % (For unstructured grid)
%% Set simulation parameters
T      = 120*second();
dT     = T/120;
dTplot = 10:20:T+10;


%% Generate grid

fault = faultLine{1};
if strcmp(typeOfGrid, 'composite')
    % Create semi-structured grid
    G = compositeGridPEBI(gridSize, [xmax, ymax], 'wellLines', wellLine, ...
                         'wellGridFactor', wellGridSize, ...
                         'mlqtMaxLevel', 2,...% 'mlqtLevelSteps', mlqtSizes,...
                         'faultLines', faultLine, 'padding', 1);
elseif strcmp(typeOfGrid, 'distmesh')
    %Create fully unstructured grid
    G = compositeGridPEBIdistmesh(gridSize, [xmax, ymax], 'wellLines', wellLine, ...
                                 'wellGridFactor', wellGridSize, ...
                                 'wellRefDist' ,wellEps, 'faultLines', faultLine, 'padding',0);
elseif strcmp(typeOfGrid, 'coarseCart')
    nx = ceil(xmax/gridSize);
    ny = ceil(ymax/gridSize);
    G = cartGrid([nx,ny],[xmax,ymax]);
    G = computeGeometry(G);
    
elseif strcmp(typeOfGrid, 'fineCart')
    gridSize = gridSize/5;
    nx = ceil(xmax/gridSize);
    ny = ceil(ymax/gridSize);
    G = cartGrid([nx,ny],[xmax,ymax]);
    G = computeGeometry(G);
end    

if strcmp(typeOfGrid,'coarseCart') || strcmp(typeOfGrid, 'fineCart')
    % Find wells OBS!!! THIS IS VERY BUGGY!!
    w1 = wellLine{1};
    w2 = wellLine{2};
    w = [w1;w2];
    D = pdist2(G.cells.centroids, w);
    [~, I] = min(D, [], 1);
    G.cells.isWell = false(G.cells.num,1);
    G.cells.isWell(I') = true(size(I'));
    % Find faults
    n = nx*100;

    fault(:,2) = fault(:,2);    
    dx = fault(2,1) - fault(1,1);
    dy = fault(2,2) - fault(1,2);
    
    vx = fault(2,:) - fault(1,:);
    spacing = linspace(0,1,n)';
    liney = fault(1,2) +ceil((dy*spacing- mod(dy*spacing, 0.5*gridSize))/gridSize)*gridSize;
    linex = fault(1,1) + dx*spacing- mod(dx*spacing, gridSize);
    line = [linex,liney];
    [line, ~, IC] = uniquetol(line,gridSize*1e-6, 'ByRows', true);
    line = line(IC,:);
    line = unique(line,'rows','stable');
    line = 0.5*(line(1:end-1,:)+line(2:end,:));
    
    D = pdist2(G.faces.centroids, line);
    [~, I] = min(D, [], 1);
    G.faces.isFault = false(G.faces.num,1);
    G.faces.isFault(I') = true(size(I'));
 
end

%Set internal boundary
Gn = makeInternalBoundary(G, find(G.faces.isFault));

% Create 3d grid
G3D = makeLayeredGrid(Gn, 5);
G3D = computeGeometry(G3D);


%% Set fluid properties
gravity reset off 

fluid = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
                        'rho', [1000, 700]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%% Set rock
rock.poro = ones(G3D.cells.num,1)*0.15;
rock.perm = ones([G3D.cells.num,1])*100*milli*darcy;


%% Create Wells
wells = find(G3D.cells.isWell);
pv = sum(poreVolume(G3D, rock));
W = verticalWell([], G3D, rock, wells(1),[], 'radius', wellGridSize/10,...
                'type', 'rate', 'val', 0.01*pv, 'Comp_i', [1,0], 'name','$I$');
W = verticalWell(W, G3D, rock, wells(2),[], 'radius', wellGridSize/10,...
                'type','bhp', 'val', 300*darcy(),'Comp_i', [0,1], 'name','$P$' );  


%% Reservoir states
state = initState(G3D, W, 0, [0.0,1]);

%% Solve
S = computeMimeticIP(G3D, rock);
state = incompMimetic(state, G3D, S, fluid, 'wells', W);

% Prepare plotting of saturations
clf;
hold on
plotGrid(G3D, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
plotWell(G3D, W, 'height', 10, 'color', 'c');
axis off equal, view([-120,30]), colormap(flipud(jet))

colorbar; hs = []; ha=[]; zoom(1.3);
 
tsave = [];
Wsave = [];
% Start the main loop
t  = 0;  plotNo = 1;
while t < T,
   state = implicitTransport(state, G3D, dT, rock, fluid, 'wells', W);

   % Check for inconsistent saturations
   assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);

   % Update solution of pressure equation.
   state = incompMimetic(state, G3D, S, fluid, 'wells', W);

    % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   tsave = [tsave; t];
   Wsave = [Wsave; state.s(wells(2:end),:)];
   
   if ( t + dT <= dTplot(plotNo)), continue, end

%    % Calculate streamlines
%    seed = (nx:nx-1:nx*ny).';
%    Sf = pollock(G3D, state, seed, 'substeps',1);
%    Sb = pollock(G3D, state, seed, 'substeps', 1,'reverse',true);
%    hf = streamline(Sf);
%    hb = streamline(Sb);
%   set([hf; hb],'color', 'k');
   % Plot saturation
   delete([hs, ha])
   hs = plotCellData(G3D, state.s(:,1), find(state.s(:,1) >= 0.0));
   ha = annotation('textbox', [0.1714 0.8214 0.5000 0.1000], 'LineStyle', 'none', ...
                   'String', ['Water saturation at ', ...
                              num2str(convertTo(t,second)), ' s']);

   plot3(fault(:,1),fault(:,2),-0.01*ones(size(fault,1),1), 'm', 'linewidth', 2);

   fig = gcf();
   set(findall(fig,'-property','FontSize'),'FontSize',14) 
   view(-120, 40), drawnow, caxis([0 1])
%    if strcmp(fileFormat, 'pdf')
%         name = strcat(typeOfGrid, num2str(convertTo(t,second)), '.pdf');
%         print('-painters', '-dpdf', '-r300', name)
%    elseif strcmp(fileFormat, 'pdfNoVector')
%         name = strcat(typeOfGrid, num2str(convertTo(t,second)), '.pdf');
%         print('-opengl', '-dpdf', '-r600', name)
%    elseif strcmp(fileFormat, 'eps')
%        name = strcat(typeOfGrid, num2str(convertTo(t,second)), '.eps');
%        print('-painters', '-depsc', '-r300', name)
%    else
%        warning('Did not recognize file type. Does not save figure')
%        plotNo = plotNo+1;
%        continue
%    end
   
   plotNo = plotNo+1;
end

close all
save(typeOfGrid)


