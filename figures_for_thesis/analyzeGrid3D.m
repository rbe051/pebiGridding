close all; clear


gridType = {'fineCartcase1','coarseCartcase1','distmeshcase1'};
legendName = {'Fine Cartesian', 'Coarse Cartesian', ...
              'Algorithm 2'};
time = [];
sat = [];
flux = [];
for i = 1:numel(gridType)
    load(gridType{i});
    time = [time,tsave];
    sat = [sat, Wsave(2:end)];
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
       plot(time(:,i), sat(:,i)*flux(i))
    end
    legend(legendName)
    %axis([0, 120, 0, 3*1.1])
    xlabel('Time (s)')
    ylabel('Oil Flux ({m^3}s^{-1})')

    fig = gcf();
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 12);
end
