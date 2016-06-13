function [state] = simulateFracturedFlowNoHybrid(G, W,rock, fluid,  endT,dt, bc)

% Transmisibility
T = computeTrans(G,rock);


%% Solve a single-phase problem, compute time of flight, and plot
stateInit = initState(G,W,0,[0 1]);
stateInit = incompTPFA(stateInit,G,T,fluid,'wells',W,'c2cTrans',T,'bc',bc);
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
if numel(W)>0
state = repmat(...
        struct('pressure',zeros(G.cells.num,1), ...
               'flux',    zeros(G.faces.num,1), ...
               's',       zeros(G.cells.num,1), ...
               'wellSol', repmat(struct('flux',[],'pressure',[]),1,numel(W)),...
               'facePressure',zeros(G.faces.num,1)),1,numSteps);
else
  state = repmat(...
        struct('pressure',zeros(G.cells.num,1), ...
               'flux',    zeros(G.faces.num,1), ...
               's',       zeros(G.cells.num,1), ...
               'facePressure',zeros(G.faces.num,1)),1,numSteps);
end
  state(1) = stateInit;
for i = 2:numSteps
    state(i) = explicitTransport(state(i-1),G,t(i),rock,fluid,'wells',W,'bc',bc);
    %state = explicitTransport_DFM(state,G,t + dt,rock,fluid,'wells',W);


    clf
    plotCellData(G,state(i).s(:,1));
    plotGrid(G,'facecolor','none')
    axis equal, axis off
    title(['Water saturation at t = ' num2str(t(i)) 's']);
    colorbar

    pause(.1)
end

end