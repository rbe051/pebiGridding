clear; close all;
addpath ../../voronoi2D/
addpath ../../../distmesh/
set(0,'defaulttextinterpreter','none')
%%
   
nF = 11;
minFracSize = 0.2;
w = cell(1,nF);
i= 1;
rng(2);
while i<=nF -1
  frac = 1.2*rand(2,2)-0.1;
  frac = min(frac,1);
  frac = max(frac,0);
  if sum(frac(1,:).^2 - frac(2,:).^2)<minFracSize^2
    w{i} = frac;
    i = i+1;
  end
end

w{6} = [0.6,0.5;1,0.1];
w{8}(2,:)= [0.2,0.8];
w{9}(2,1)=0.4;
w{nF} = [0.6,0.6;0.95,1];
gs = 0.2;
wGf = 1/8;
rng(1);

ap = 1e-3;
resDim = 1000;
d = {@(p) ap/resDim + 0*ap/100/resDim*rand(size(p,1),1)};


close all
figure()

plotLinePath(w)
%% Create Grid
[G,Gn] = pebiGridForCompareFractures(gs,[1,1], 'wellLines',w,...
                      'protLayer',true,'protD',d,'wellRefinement',true,...
                      'wellGridFactor',wGf,'epsilon',0.15);

G.nodes.coords = G.nodes.coords*resDim;
Gn.nodes.coords = Gn.nodes.coords*resDim;
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
aperture = ap;%dist(Gn.faces.centroids(Gn.faces.tag,:));
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
rockn.perm = 1/100*milli * darcy * ones(nCells,2); % Edge centered
rockn.poro = 0.001 * ones(nCells,1);               % Edge centered
rock.perm = 1/100*milli * darcy * ones(cells,2);   % Centroid centered
rock.poro =  0.001 * ones(cells,1);                % Centroid centered


% The fracture permeability is computed from the parallel plate assumption,
% which states that the permeability is aperture squared divided by 12.
poro = 0.5;
rockn.perm(hybridInd,:) = repmat(aperture.^2/12,sum(hybridInd>0),2);    % Edge centered
rockn.poro(hybridInd) = poro;                                            % Edge centered
rock.perm(G.cells.tag,:) = repmat(aperture.^2/12,sum(G.cells.tag),2);   % Centroid centered
rock.poro(G.cells.tag) = poro;                                           % Centroid centered

% Fix porosity for intersection cells
medVol = median(G.cells.volumes(G.cells.tag));
toBig = G.cells.volumes(G.cells.tag)>10*medVol;
c = find(G.cells.tag);
rock.poro(c(toBig)) = poro*medVol./G.cells.volumes(c(toBig));


% Create fluid object. The fluids have equal properties, so this ammounts
% to the injection of a tracer
fluid = initSimpleFluid('mu' , [   1,  1]*centi*poise     , ...
    'rho', [1014, 859]*kilogram/meter^3, ...
    'n'  , [   1,   1]);



%% Add boundary condition
% BC Edge centered
tol = 1e-6;
lowerBn = Gn.faces.centroids(:,2)<tol;
upperBn = Gn.faces.centroids(:,2)>resDim - tol;

bcVal = 0.1/resDim/day();
area = G.faces.areas;
arean = Gn.faces.areas;

bcn=addBC([],find(lowerBn),'flux',bcVal*arean(lowerBn));
bcn=addBC(bcn,find(upperBn),'pressure',0);
bcn.sat = zeros(size(bcn.face,1),2);
bcn.sat(:,1)=1;
% BC centroid  centered
lowerB = G.faces.centroids(:,2)<tol;
upperB = G.faces.centroids(:,2)>resDim - tol;
bc     =addBC([],find(lowerB),'flux',bcVal*area(lowerB));
bc     =addBC(bc,find(upperB),'pressure',0);
bc.sat = zeros(size(bc.face,1),2);
bc.sat(:,1)=1;
W = {};
Wn = {};
%% Solve
endTime = 0.2*sum(poreVolume(G,rock) / (bcVal*resDim));
numSteps = 100;
dt = endTime / numSteps;
  
staten = simulateFracturedFlow(Gn,Wn,rockn, fluid, endTime, dt,bcn);
state = simulateFracturedFlowNoHybrid(G,W,rock, fluid, endTime, dt,bc);

%% Save workspace
save('fracturedReservoir')
%% Plot end sol
%close all
t = 0:dt:endTime;
for i = 1:numel(staten)
  figure(1)
  plotCellData_DFM(Gn,staten(i).s(:,1));
  plotFractures(Gn,hybridInd,staten(i).s(:,1));
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
  
  pause(0.01)
  clf
  figure(1)
  clf
end



%% Plot flow over upper boundary
t = 0:dt:endTime;
flow = zeros(numel(t),1);
flown = zeros(numel(t),1);
c = sum(G.faces.neighbors(upperB,:),2); %should give cell number since other neighbor is 0
cn = sum(Gn.faces.neighbors(upperBn,:),2);
for i = 1:numel(staten)
  sn    = staten(i).s(:,2);
  fluxn = staten(i).flux;
  s     =  state(i).s(:,2);
  flux  =  state(i).flux;
  
  
  flow(i)  = sum(s(c).*flux(upperB).*dt);
  flown(i) = sum(sn(cn).*fluxn(upperBn).*dt);
end
figure(); hold on
plot(t/day(),cumsum(flow))
plot(t/day(),cumsum(flown))


%% Plot error
% find reservoir cells
close all
c = ~(G.cells.tag|G.cells.protectionCells);
intF = Gn.faces.tags>0;
intC = Gn.faces.neighbors(intF,:);
cn = true(Gn.cells.num,1);
cn(intC) = false;
cn(Gn.cells.hybrid>0) = false;

err = zeros(sum(c),numel(state));
for i = 1:numel(staten)
  p = state(i).pressure;
  pn = staten(i).pressure;
  
  s     =  state(i).s(:,1);
  sn    = staten(i).s(:,1);
    
  err(:,i) = (p(c) - pn(cn))/max(p(c));
  errs(:,i) = (s(c) - sn(cn));
end

plotCellData(G,err(:,1),c)

figure(1)
plot(t/day(), sqrt(sum(errs.^2,1)))
%figure()
%for i =1:numel(t)
%  plotCellData(G,errs(:,i),c)
%  colorbar()
%  pause(0.1)
%end
figure(2)
plotAtDay = 500;
[~, id] = min(abs(t/day -plotAtDay));
plotCellData(G,errs(:,id),c)
colorbar
axis equal off tight

figure(3);
plotCellData_DFM(Gn,staten(id).s(:,1))
plotGrid_DFM(Gn,'facecolor','none','edgecolor',[0.1,.1,.1])
colorbar
axis equal off tight

%% Save
%figure(1)
%print('../../../../master/thesis/fig/ch06/inkscape/diffInSatErrorPlot','-dsvg')
figure(2)
print('../../../../master/thesis/fig/ch06/inkscape/diffInSatGridPlot','-dsvg')
figure(3)
print('../../../../master/thesis/fig/ch06/inkscape/saturationFracturedRes','-dsvg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot fractures
close all
figure(1);hold on
color = get(gca,'ColorOrder');
plotLinePath(w,'linewidth',2,'color',color(4,:))
axis equal off tight
plot([0,0,1,1,0], [0,1,1,0,0],'k--')
%% Save fractures
figure(1)
print('../../../../master/thesis/fig/ch06/fracturedRes','-depsc')




