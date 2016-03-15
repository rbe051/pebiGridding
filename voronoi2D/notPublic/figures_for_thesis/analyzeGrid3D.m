close all; clear


gridType = {'fineCart','coarseCart','composite', 'distmesh'};
legendName = {'Fine Cartesian', 'Coarse Cartesian', 'Algorithm 1', ...
              'Algorithm 2'};
time = [];
sat = [];
flux = [];
for i = 1:numel(gridType)
    load(gridType{i});
    time = [time,tsave];
    sat = [sat, Wsave(2:end,1)];
    newFlux = [];
    for j = 2:numel(state.wellSol)
        newFlux = [newFlux,abs(sum(state.wellSol(j).flux))];
    end    
    flux = [flux; newFlux];
    G3D.cells.num
end

    

for j = 1:2:size(sat,2)-1
    figure()
    hold on
    for i = 1:numel(gridType)
       plot(time(2:end,i), (1-sat(:,i))*flux(i))
    end
    legend(legendName)
    %axis([0, 120, 0, 3*1.1])
    xlabel('Time (s)')
    ylabel('Oil Flux ({m^3}s^{-1})')

    fig = gcf();
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 12);
end

axis([0,120,0,3.5])