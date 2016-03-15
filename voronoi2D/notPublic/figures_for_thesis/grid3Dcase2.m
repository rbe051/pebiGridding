%% Grid with single source centered at 0.5,0.5
close all; clear;

addpath('../', '../mrstTweak/')
mrstModule add mimetic
mrstModule add distmesh


typeOfGrid = 'composite';
name = 'case2';
name = strcat(name,typeOfGrid);

fileFormat = 'pdfNoVector';

%% grid parameters


%[nx,ny,nz] = deal(160,160,10); %fine cart 
%[nx,ny,nz] = deal(10,10,10);%coarse cart
[nx,ny,nz] = deal(20,20,10);%distmesh
%[nx,ny,nz] = deal(10,10,10);%composite
[Lx,Ly,Lz] = deal(200,200,50);

gridSize = norm([Lx,Ly])/norm([nx,ny]);                 % Set size of grid cells

faultLine = {[40,60;160,180], [40,120;160,40]};
wellLine = {[100,50]};                   % Set source center

mlqtMax = 0;                        % Set number of reminement levels
wellGridSize = 0.5/2^mlqtMax;       % Set relative well grid size
mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                    % Set distance around wells to be
                                    % refined.

wellEps = (norm(Lx,Ly))*1/10;           % Size around wells to be refined
                                       % (For unstructured grid)

%% Generate grid


if strcmp(typeOfGrid, 'composite')
    % Create semi-structured grid
    G = compositeGridPEBI(gridSize, [Lx, Ly], 'wellLines', wellLine, ...
                         'wellGridFactor', wellGridSize, ...
                         'mlqtMaxLevel', mlqtMax,...% 'mlqtLevelSteps', mlqtSizes,...
                         'faultLines', faultLine, 'padding', 1);
elseif strcmp(typeOfGrid, 'distmesh')
    %Create fully unstructured grid
    G = compositeGridPEBIdistmesh(gridSize, [Lx, Ly], 'wellLines', wellLine, ...
                                 'wellGridFactor', wellGridSize, ...
                                 'wellRefDist' ,wellEps, 'faultLines', faultLine, 'padding',0);
else
    G = cartGrid([nx,ny] , [Lx,Ly]);

end    
    G = computeGeometry(G);
if strcmp(typeOfGrid,'coarseCart') || strcmp(typeOfGrid, 'fineCart')
    % Find wells OBS!!! THIS IS VERY BUGGY!!
    w = wellLine{1};

    D = pdist2(G.cells.centroids, w);
    [~, I] = min(D, [], 1);
    G.cells.isWell = false(G.cells.num,1);
    G.cells.isWell(I') = true(size(I'));
    
    % Find faults
    n = nx*100;
    G.faces.isFault = false(G.faces.num,1);
    for i = 1:numel(faultLine)
        fault = faultLine{i};
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
        G.faces.isFault(I') = true(size(I'));
    end
end
if strcmp(typeOfGrid,'composite') || strcmp(typeOfGrid, 'distmesh')
    
    % Find wells OBS!!! THIS IS VERY BUGGY!!
    w = wellLine{1};

    D = pdist2(G.cells.centroids, w);
    [~, I] = min(D, [], 1);
    G.cells.isWell = false(G.cells.num,1);
    G.cells.isWell(I') = true(size(I'));
    % Find wells OBS!!! THIS IS VERY BUGGY!!
    % Find faults
    n = nx;
    for i = 2:numel(faultLine)
        fault = faultLine{i};
        fault(:,2) = fault(:,2);    
        dx = fault(2,1) - fault(1,1);
        dy = fault(2,2) - fault(1,2);
        spacing = linspace(0.28,0.35,n)';
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
G = makeLayeredGrid(Gn, nz);
G.nodes.coords(:,3) = G.nodes.coords(:,3)*Lz/nz;
G = computeGeometry(G);

%% Set rock
rock.poro = ones(G.cells.num,1)*0.3;
rock.perm = ones([G.cells.num,1])*30*milli*darcy;

%% Porve volume vs pressure
cr = 1e-6/barsa;
p_r = 200*barsa;
pv_r = poreVolume(G,rock);

pv = @(p) pv_r.*exp(cr* (p - p_r));

%% Density at reference pressure
mu    = 5*centi*poise;
c     = 1e-3/barsa;
rho_r = 850*kilogram/meter^3;
rhoS  = 750*kilogram/meter^3;
rho   = @(p) rho_r.*exp(c*(p-p_r));


%% Create Wells
wells = find(G.cells.isWell);
W = verticalWell([], G, rock, wells(1),[], 'name','$P$');
%% Find hydrestatic preassure
gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil = ode23(@(z,p) g.*rho(p), [z_0, z_max],p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)),[],1);

%% Discretization
N = double(G.faces.neighbors);
intInx = all(N ~= 0,2);
N = N(intInx,:);

n = size(N,1);
C= sparse([(1:n)';  (1:n)'],N, ones(n,1)*[-1,1], n, G.cells.num);
grad = @(x) C*x;
div = @(x) -C'*x;
avg = @(x) 0.5*(x(N(:,1)) + x(N(:,2)));

%% Transmitibilities
hT = computeTrans(G,rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T = 1./accumarray(cf, 1./hT, [nf,1]);
T = T(intInx);

%% Residual eequations
gradz = grad(G.cells.centroids(:,3));
v = @(p) -(T/mu).*(grad(p)-g*avg(rho(p)).*gradz);
presEq=@(p,p0,dt)  (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0)) ...
                   + div(avg(rho(p)).*v(p) );




%% well model
wc = W(1).cells;
WI = W(1).WI;
dz = W(1).dZ;

p_conn = @(bhp) bhp + g*dz.*rho(bhp);
q_conn = @(p,bhp) WI.*(rho(p(wc))/mu).*(p_conn(bhp) - p(wc));
rateEq = @(p, bhp, qS) qS - sum(q_conn(p, bhp))/rhoS;

ctrlEq = @(bhp) bhp - 100*barsa;

%% Solve
[p_ad,bhp_ad,qS_ad] =initVariablesADI(p_init,p_init(wc(1)), 0);

nc = G.cells.num;
[pIx, bhpIx, qSIx] = deal(1:nc, nc+1,nc+2);
numSteps = 52;
totTime = 2*year;
dt = totTime/numSteps;
tol = 1e-5;
maxits = 10;
plotTimes = (50*day:50*day:totTime);

sol  = repmat(struct('time', [], 'pressure', [], 'bhp', [], ...
                    'qS', [], 'rate', []), [numSteps + 1, 1]);
sol(1) =      struct('time',0, 'pressure', double(p_ad), ...
                     'bhp', double(bhp_ad), 'qS', double(qS_ad),...
                     'rate', double(rateEq(p_ad,bhp_ad,qS_ad)));

%% INitial plot
show = true(G.cells.num,1);
% cellInx = sub2ind(G.cartDims, ...
%     [I-1;I-1;I;I  ;I(1:2) -1], ...
%     [J  ;J  ;J;J  ;nperf + [2;2]],...
%     [K-1;K  ;K;K-1;K(1:2)-[0;1]]);
%show(cellINx) = false;


% Prepare plotting of pressure
clf;
hold on
plotCellData(G, double(p_ad)/barsa);
plotWell(G, W, 'height', 50, 'color', 'c');
for i = 1:numel(faultLine)
  line = faultLine{i};
   plot3(line(:,1),line(:,2),-0.01*ones(size(line,1),1), 'm', 'linewidth', 2);
end

axis off equal, view([-120,30]), colormap(flipud(jet))
colorbar; hs = []; ha=[]; zoom(1.3);
pause()
%% The time loop
t = 0; step = 0; plotNo=1;
while t < totTime
    t = t + dt; step = step + 1;
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n',...
            step, convertTo(t-dt,day), convertTo(t,day));
    % Newton loop
    resNorm = 1e99;
    p0 = double(p_ad);
    nit = 0;
    while (resNorm > tol) && (nit <=maxits)
        eq1 = presEq(p_ad, p0, dt);
        eq1(wc) = eq1(wc) -q_conn(p_ad, bhp_ad);
        eqs = {eq1, rateEq(p_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};
        eq = cat(eqs{:});
        
        J = eq.jac{1};
        res = eq.val;
        upd = -(J\res);
        p_ad.val = p_ad.val + upd(pIx);
        bhp_ad.val = bhp_ad.val + upd(bhpIx);
        qS_ad.val = qS_ad.val +upd(qSIx);
        
        
        resNorm = norm(res);
        nit = nit + 1;
        fprintf(' Iteration %3d: Res = %.4e\n', nit, resNorm);
        
    end

    if nit >maxits, error('Newton solves did not converge')
    else 
        sol(step+1) = struct('time', t,'pressure', double(p_ad), ...
                             'bhp', double(bhp_ad),'qS', double(qS_ad),...
                             'rate', double(rateEq(p_ad,bhp_ad,qS_ad*0)));
    end
    
    %% plotting
    if (plotNo<=numel(plotTimes) &&  t <= plotTimes(plotNo)), continue, end

%    % Calculate streamlines
%    seed = (nx:nx-1:nx*ny).';
%    Sf = pollock(G3D, state, seed, 'substeps',1);
%    Sb = pollock(G3D, state, seed, 'substeps', 1,'reverse',true);
%    hf = streamline(Sf);
%    hb = streamline(Sb);
%   set([hf; hb],'color', 'k');
   % Plot saturation
   delete([hs, ha])
   hs = plotCellData(G, double(p_ad)/barsa);
   ha = annotation('textbox', [0.1714 0.8214 0.5000 0.1000], 'LineStyle', 'none', ...
                   'String', ['Pressure at ', ...
                              num2str(round(convertTo(t,day))), ' days']);
   fig = gcf();
   set(findall(fig,'-property','FontSize'),'FontSize',14) 
   view(-120, 30), drawnow, caxis([100 200])
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

save(name,'sol','rateEq');

