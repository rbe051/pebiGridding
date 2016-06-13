function [state] = simulateFracturedFlow(G, W,rock, fluid,  endT,dt,bc)
T = computeTrans_DFM(G,rock,'hybrid',true);

% Transmissibilites for fracture-fracture connections are computed in a
% separate file
[G,T2] = computeHybridTrans(G,T);
hybridInd = find(G.cells.hybrid);


%% Solve a single-phase problem, compute time of flight, and plot
stateInit = initState(G,W,0,[0 1]);
stateInit = incompTPFA_DFM(stateInit,G,T,fluid,'wells',W,'c2cTrans',T2,'bc',bc);
%state = incompTPFA_DFM(state,G,T,fluid,'wells',W,'c2cTrans',T2);



% Then solve tracer transport

t = 0;

% End of simulation


% Since the two fluids have equal properties, the pressure solution is time
% independent, and the transport equation can be solved for the entire
% simulation time at once. For visualization purposes, we split the
% interval anyhow




t = 0:dt:endT;
numSteps = numel(t);
if numel(W) >0
state = repmat(...
        struct('pressure',zeros(G.cells.num,1), ...
               'flux',    zeros(G.faces.num,1), ...
               's',       zeros(G.cells.num,1), ...
               'wellSol', repmat(struct('flux',[],'pressure',[]),1,numel(W)),...
               'fluxc2c', zeros(size(G.cells.neighbors,1),1),...
               'facePressure',zeros(G.faces.num,1)),1,numSteps);
else
state = repmat(...
        struct('pressure',zeros(G.cells.num,1), ...
               'flux',    zeros(G.faces.num,1), ...
               's',       zeros(G.cells.num,1), ...
               'fluxc2c', zeros(size(G.cells.neighbors,1),1),...
               'facePressure',zeros(G.faces.num,1)),1,numSteps);
end  
state(1) = stateInit;
for i = 2:numSteps
    state(i) = explicitTransport_DFM(state(i-1),G,t(i),rock,fluid,'wells',W,'bc',bc);
    %state = explicitTransport_DFM(state,G,t + dt,rock,fluid,'wells',W);

    clf
    plotCellData_DFM(G,state(i).s(:,1));
    plotFractures(G,hybridInd,state(i).s(:,1));
    plotGrid_DFM(G,'facecolor','none')
    axis equal, axis off
    title(['Water saturation at t = ' num2str(t(i)) 's']);
    colorbar

    pause(.1)
end

end