clear; close all;
addpath ../../voronoi2D/
addpath ../../../distmesh/
set(0,'defaulttextinterpreter','none')
%%
xmax = 100;
x = 1:0.01:10-1;
y = 2.5+4*sin(pi*x/10);
w = {[x',y']};

gs = 2;
wGf = 1/8;
rng(1);
d = {@(p) 1e-2+ 0.05*rand(size(p,1),1)};


%% Create Grid
[G,Gn] = pebiGridForCompareFractures(gs,[10,10], 'wellLines',w,...
                      'protLayer',true,'protD',d,'wellRefinement',true,...
                      'wellGridFactor',wGf,'epsilon',2);

G.nodes.coords = G.nodes.coords*xmax/10;
Gn.nodes.coords = Gn.nodes.coords*xmax/10;
G = computeGeometry(G);
Gn = computeGeometry(Gn);
%% Plot Grid

plotGrid(G)
axis equal
figure
plotGrid(Gn)
axis equal
%%   Assign aperture
% For edge Grid:
try
   require dfm incomp
catch
   mrstModule add dfm incomp
end

rng(1);
dist = d{1};
Gn.faces.tags = Gn.faces.tag;
apt = zeros(Gn.faces.num,1);
aperture = dist(Gn.faces.centroids(Gn.faces.tag,:))*xmax/10;
apt(Gn.faces.tags) = aperture;


%   Add the hybrid cells
Gn = addhybrid(Gn,Gn.faces.tag > 0,apt);

% Plot the grid and the fracture
figure
plotGrid_DFM(Gn)
plotFractures(Gn)
axis equal, axis off

%% Add physical parameters

% Find indices of hybrid cells
hybridInd = find(Gn.cells.hybrid);

nCells = Gn.cells.num;
cells = G.cells.num;
% Define permeability and porosity
rockn.perm = milli * darcy * ones(nCells,2); % Edge centered
rockn.poro = 0.01 * ones(nCells,1);          % Edge centered
rock.perm = milli * darcy * ones(cells,2);   % Centroid centered
rock.poro =  0.01 * ones(cells,1);           % Centroid centered


% The fracture permeability is computed from the parallel plate assumption,
% which states that the permeability is aperture squared divided by 12.


rockn.perm(hybridInd,:) = aperture.^2/12*[1,1];    % Edge centered
rockn.poro(hybridInd) = 0.5;                 % Edge centered
rock.perm(G.cells.tag,:) = aperture.^2/12*[1,1];    % Centroid centered
rock.poro(G.cells.tag) = 0.5;               % Centroid centered
% Create fluid object. The fluids have equal properties, so this ammounts
% to the injection of a tracer
fluid = initSimpleFluid('mu' , [   1,  1]*centi*poise     , ...
    'rho', [1014, 859]*kilogram/meter^3, ...
    'n'  , [   1,   1]);



%% Add boundary condition
bcType = 'flux';
bcVal = 0;
%% Add wells
wellRate = 1;
bhp = 1e8;
% Injection well in lower left corner, production in upper right

w1 = 1;
w2 = sum(Gn.faces.tag);

Wn = addWell([],Gn,rockn,w1,'type','bhp','val',bhp,'comp_i',[1,0],'InnerProduct','ip_tpf','radius',0.01);
Wn = addWell(Wn,Gn,rockn,w2,'type','bhp','val',0,'comp_i',[1 0],'InnerProduct','ip_tpf','radius',0.01);

W  = addWell([],G,rock,w1,'type','bhp','val',bhp,'comp_i',[1,0],'InnerProduct','ip_tpf','radius',0.01);
W = addWell(W,G,rock,w2,'type','bhp','val',0,'comp_i',[1 0],'InnerProduct','ip_tpf','radius',0.01);
%% Solve
endTime = 1/12*hour*sum(0.1 * poreVolume(G,rockn) / wellRate);
numSteps = 100;
dt = endTime / numSteps;
  
staten = simulateFracturedFlow(Gn,Wn,rockn, fluid, endTime, dt,bcType,bcVal);
state = simulateFracturedFlowNoHybrid(G,W,rock, fluid, endTime, dt,bcType,bcVal);

%% Plot end sol
close all
t = 0:dt:endTime;
for i = 1:numel(staten)
  figure(1)
  plotCellData_DFM(Gn,staten(i).s(:,1));
  plotFractures(Gn,hybridInd,staten(end).s(:,1));
  title(['Water saturation at t = ' num2str(t(i)) 's']);
  plotGrid_DFM(Gn,'facecolor','none')
  axis equal, axis off

  colorbar

  figure(2)
  plotCellData(G,state(i).s(:,1));
  plotGrid(G,'facecolor','none')
  title(['Water saturation at t = ' num2str(t(i)) 's']);
  axis equal, axis off
  colorbar
  
  pause()
  clf
  figure(1)
  clf
end


%% Plot Production
fluxn = zeros(numel(staten),1);
satn = zeros(numel(staten),2);
flux = zeros(numel(state),1);
sat = zeros(numel(state),2);
t = 0:dt:endTime;
for i = 1:numel(staten)
  fluxn(i) = staten(i).wellSol(2).flux;
  satn(i,:) = staten(i).s(Wn(2).cells,:);
  flux(i) = state(i).wellSol(2).flux;
  sat(i,:) = state(i).s(W(2).cells,:);
end
figure(); hold on
color = get(gca,'ColorOrder');
plot(t/hour(), satn(:,1),'color',color(4,:),'linewidth',1.2)
plot(t/hour(), sat(:,1),'color',color(5,:),'linewidth',1.2)

xlabel('Number of Hours')
ylabel('Water Saturation')
%legend('Edge Centered Fracture','Centroid Centered Fracture','Location','northwest')
%% Save 
print('../../../../master/thesis/fig/ch06/compareFracGridMethodSingleFault','-dsvg')


%% Plot reseroir at given times
close all
pltAtHour = [3,6];
name = {};
for i = 1:numel(pltAtHour)
  j = round(pltAtHour(i)*hour/dt);
  figure(i)
  plotCellData_DFM(Gn,staten(j).s(:,1));
  plotFractures(Gn,hybridInd,staten(end).s(:,1));
  name{i} = (['EdgeWaterSaturationAt' num2str(pltAtHour(i)) 'h']);
  plotGrid_DFM(Gn,'facecolor','none')
  axis equal, axis off
  colorbar
  
  figure(i+numel(pltAtHour))
  plotCellData(G,state(j).s(:,1));
  plotGrid(G,'facecolor','none')
  name{i+numel(pltAtHour)} = (['CentroidWaterSaturationAt' num2str(pltAtHour(i)) 'h']);
  axis equal, axis off
  colorbar
end

%% Save
for i = numel(name)
  figure(i)
  fName = strcat('../../../../master/thesis/fig/ch06/inkscape/singleFrac', name{i});
  print(fName,'-dsvg')
end






