clear;close all
%addpath ../../voronoi2D/
%% Create fault and well

w = {[0.1,0.3;0.9,0.3], ...
     [0.1,0.7;0.9,0.7],};

f = {[0.3,0.05;0.2,0.5],...
     [0.35,0.1;0.45,0.65],...
     [0.4,0.1;0.63,0.8],...
     [0.55,0.05;0.48,0.6],...
     [0.65,0.1;0.585,0.59],...
     [0.7,0.05;0.7,0.5],...
     [0.8,0.05;0.75,0.4], ...
     [0.3,0.5;0.2,0.9],...
     [0.4,0.5;0.4,0.95],...
     [0.3,0.45;0.6,0.95],...
     [0.55,0.4;0.5,0.9],...
     [0.57,0.4;0.7,0.9],...
     [0.7,0.5;0.7,0.9],...
     [0.8,0.4;0.75,0.8]};

rng(8);
%f = cellfun(@(c) c + 0.1*rand(size(c)),f,'un',false);
%% Plot
figure()
plotLinePath(f)
plotLinePath(w)
%% create grid
Nx = 1; Ny = 1;
gSize = 10;
rGs = 0.2*Nx;
wGf = 1/15;
fGf = 1/1.5;
eps = 1/4;
G = pebiGrid(rGs, [Nx,Ny],'wellLines',w,'faultLines',f, ...
                      'wellGridFactor',wGf, 'faultGridFactor',fGf,...
                      'wellRefinement',true,'epsilon',eps,'circlefactor',0.57);

G.nodes.coords = G.nodes.coords*gSize;
G = computeGeometry(G);
%% Solve
try
   require dfm incomp
catch
   mrstModule add dfm incomp
end


%%   Assign aperture
G.faces.tags = G.faces.tag;
apt = zeros(G.faces.num,1);
aperture = 0.001;
apt(G.faces.tags) = aperture;


%   Add the hybrid cells
G = addhybrid(G,G.faces.tag > 0,apt);

% Plot the grid and the fracture
figure
plotGrid_DFM(G)
plotFractures(G)
axis equal, axis off

%% Add physical parameters

% Find indices of hybrid cells
hybridInd = find(G.cells.hybrid);

nCells = G.cells.num;

% Define permeability and porosity
rock.perm = milli * darcy * ones(nCells,2);
rock.poro = 0.01 * ones(nCells,1);

% The fracture permeability is computed from the parallel plate assumption,
% which states that the permeability is aperture squared divided by 12.
rock.perm(hybridInd,:) = aperture^2/12;
rock.poro(hybridInd) = 0.5;
rock.perm(G.cells.tag) = 10*aperture^2/12;
rock.poro(G.cells.tag) = 1;
% Create fluid object. The fluids have equal properties, so this ammounts
% to the injection of a tracer
fluid = initSimpleFluid('mu' , [   1,  1]*centi*poise     , ...
    'rho', [1014, 859]*kilogram/meter^3, ...
    'n'  , [   1,   1]);
% Compute TPFA transmissibilities
T = computeTrans_DFM(G,rock,'hybrid',true);

% Transmissibilites for fracture-fracture connections are computed in a
% separate file
[G,T2] = computeHybridTrans(G,T);


%% Add boundary condition
tol = 1e-6;
f = any(G.faces.centroids<=tol,2);
f = f | any(G.faces.centroids>=gSize-tol,2);
bc=addBC([],find(f),'flux',0);
bc.sat = zeros(size(bc.face,1),2);
bc.sat(:,1)=1;
%% Add wells
wellRate = 1;
bhp = 1e6;
% Injection well in lower left corner, production in upper right
w = find(G.cells.tag);
wCent = G.cells.centroids(w,:);
dist =  sum(wCent.^2,2);
[~, Imin]= min(dist);
[~, Imax] = max(dist);
w1 = w(Imin);
w2 = w(Imax);

W = addWell([],G,rock,w1,'type','bhp','val',bhp,'comp_i',[1,0],'InnerProduct','ip_tpf','radius',0.01);
W = addWell(W,G,rock,w2,'type','bhp','val',-bhp,'comp_i',[1 0],'InnerProduct','ip_tpf','radius',0.01);

%% Solve a single-phase problem, compute time of flight, and plot
state = initState(G,W,0,[0 1]);
state = incompTPFA_DFM(state,G,T,fluid,'wells',W,'c2cTrans',T2,'bc',bc);
%state = incompTPFA_DFM(state,G,T,fluid,'wells',W,'c2cTrans',T2);

figure
plotCellData_DFM(G,state.pressure)
plotFractures(G,hybridInd,state.pressure)
axis equal, axis off
colorbar

%% Then solve tracer transport

t = 0;

% End of simulation
endTime = 100*sum(0.1 * poreVolume(G,rock) / wellRate);

% Since the two fluids have equal properties, the pressure solution is time
% independent, and the transport equation can be solved for the entire
% simulation time at once. For visualization purposes, we split the
% interval anyhow
numSteps = 5;
dt = endTime / numSteps;

iter = 1;

while t < endTime
    state = explicitTransport_DFM(state,G,t + dt,rock,fluid,'wells',W,'bc',bc);
    %state = explicitTransport_DFM(state,G,t + dt,rock,fluid,'wells',W);

    iter = iter + 1;
    t = t + dt;

    clf
    plotCellData_DFM(G,state.s(:,1));
    plotFractures(G,hybridInd,state.s(:,1));
    plotGrid_DFM(G,'facecolor','none')
    axis equal, axis off
    title(['Water saturation at t = ' num2str(t) 's']);
    colorbar

    pause(.1)
end
