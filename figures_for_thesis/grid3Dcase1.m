%% Grid with single source centered at 0.5,0.5
close all; clear;
addpath('../', '../mrstTweak/')
mrstModule add mimetic
mrstModule add distmesh

typeOfGrid = 'coarseCart';
name = strcat(typeOfGrid, 'case1');
fileFormat = 'pdfNoVector';

%% grid parameters
%[nx,ny,nz] = deal(120,60,10); %fine cart 
[nx,ny,nz] = deal(20,10,10);%coarse cart
%[nx,ny,nz] = deal(10,0,10);%distmesh
%[nx,ny,nz] = deal(10,10,10);%composite
[Lx,Ly,Lz] = deal(400,200,50);

gridSize = norm([Lx,Ly])/norm([nx,ny]);

faultLine = {[40,60;360,180], [40,160;250,40]};
wellLine = {[40,100], [350, 150], [160,140]};                   % Set source center

mlqtMax = 2;                        % Set number of reminement levels
wellGridSize = 0.8/2^mlqtMax;       % Set relative well grid size
mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                    % Set distance around wells to be
                                    % refined.

wellEps = sqrt(gridSize)*8;      % Size around wells to be refined8
                                 % (For unstructured grid)
%% Set simulation parameters
T      = 12*hour;
dT     = T/120;
dTplot = 10:20:T+10;


%% Generate grid


if strcmp(typeOfGrid, 'composite')
    % Create semi-structured grid
    G = compositeGridPEBI(gridSize, [Lx, Ly], 'wellLines', wellLine, ...
                         'wellGridFactor', wellGridSize, ...
                         'mlqtMaxLevel', 2,...% 'mlqtLevelSteps', mlqtSizes,...
                         'faultLines', faultLine, 'padding', 1);
elseif strcmp(typeOfGrid, 'distmesh')
    %Create fully unstructured grid
    G = compositeGridPEBIdistmesh(gridSize, [Lx, Ly], 'wellLines', wellLine, ...
                                 'wellGridFactor', wellGridSize, ...
                                 'wellRefDist' ,wellEps, 'faultLines', faultLine,...
                                 'faultGridFactor',1/3,'circleFactor',0.6,'padding',0);
elseif strcmp(typeOfGrid, 'coarseCart')
    G = cartGrid([nx,ny],[Lx,Ly]);
    G = computeGeometry(G);
    
elseif strcmp(typeOfGrid, 'fineCart')
    G = cartGrid([nx,ny],[Lx,Ly]);
    G = computeGeometry(G);
end    
    G = computeGeometry(G);
if strcmp(typeOfGrid,'coarseCart') || strcmp(typeOfGrid, 'fineCart')
    % Find wells OBS!!! THIS IS VERY BUGGY!!
    w = [];
    for i = 1:numel(wellLine)
        w = [w; wellLine{i}];
    end
    D = pdist2(G.cells.centroids, w);
    [~, I] = min(D, [], 1);
    G.cells.isWell = false(G.cells.num,1);
    G.cells.isWell(I') = true(size(I'));
    
    % Find faults
    n = nx*1000;
    G.faces.isFault = false(G.faces.num,1);
    for i = 1:numel(faultLine)
        fault = faultLine{i};
        fault(:,2) = fault(:,2);    
        dx = fault(2,1) - fault(1,1);
        dy = fault(2,2) - fault(1,2);

        vx = fault(2,:) - fault(1,:);
        spacing = linspace(0,1,n)';
        sizey = Ly/ny;
        sizex = Lx/nx;
        liney = fault(1,2) +ceil((dy*spacing- mod(dy*spacing, 0.5*sizey))/sizey)*sizey;
        linex = fault(1,1) + ceil((dx*spacing- mod(dx*spacing, 0.5*sizex))/sizex)*sizex;%dx*spacing- mod(dx*spacing, gridSize);
        line = [linex,liney];
        [line, ~, IC] = uniquetol(line,gridSize*1e-6, 'ByRows', true);
        line = line(IC,:);
        line = unique(line,'rows','stable');
        line = 0.5*(line(1:end-1,:)+line(2:end,:));

        D = pdist2(G.faces.centroids, line);
        [~, I] = min(D, [], 1);
        G.faces.isFault(I') = true(size(I'));
    end
end
if strcmp(typeOfGrid,'composite') || strcmp(typeOfGrid, 'distmesh')
    % Find wells OBS!!! THIS IS VERY BUGGY!!
    % Find faults
    n = nx;
    for i = 2:numel(faultLine)
        fault = faultLine{i};
        fault(:,2) = fault(:,2);    
        dx = fault(2,1) - fault(1,1);
        dy = fault(2,2) - fault(1,2);
        spacing = linspace(0.4,0.65,n)';
        linex = fault(1,1) + dx*spacing;
        liney = fault(1,2) + dy*spacing;
        line = [linex,liney];
        [line, ~, IC] = uniquetol(line,gridSize*1e-6, 'ByRows', true);
        line = line(IC,:);
        line = unique(line,'rows','stable');
        line = 0.5*(line(1:end-1,:)+line(2:end,:));
        plotGrid(G);
        hold on
        plot(line(:,1), line(:,2),'o')
        D = pdist2(G.faces.centroids, line);
        [~, I] = min(D, [], 1);
        G.faces.isFault(I') = true(size(I'));
    end
end

% Plot grid
plotGrid(G);
hold on
plotFault(G,'color','r')
for i = 1:numel(faultLine)
  line = faultLine{i};
  plot(line(:, 1), line(:, 2),'color','m');
end
pause()

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
                'type', 'rate', 'val', 3, 'Comp_i', [1,0], 'name','$I$');
for i = 2:numel(wells)
    W = verticalWell(W, G3D, rock, wells(i),[], 'radius', wellGridSize/10,...
                    'type','bhp', 'val', 300*darcy(),'Comp_i', [0,1], 'name','$P$' );  
end

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
   wellProd = [];
   for i = 2:numel(W)
       wellCells = W(i).cells;
       wellProd = [wellProd, sum(state.s(wellCells,2).*state.wellSol(i).flux)];
   end
   
   Wsave = [Wsave; wellProd];
   
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
    
   for i = 1:numel(faultLine)
       fault = faultLine{i};
       plot3(fault(:,1),fault(:,2),-0.01*ones(size(fault,1),1), 'm', 'linewidth', 2);
   end
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
save(name)


