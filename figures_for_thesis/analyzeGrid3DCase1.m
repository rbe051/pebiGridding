close all; clear


gridType = {'fineCart','coarseCart','composite', 'distmesh'};%
legendName = { 'Fine Cartesian','Coarse Cartesian','Algorithm 1', ... %
              'Algorithm 2'};
time = [];
sat = [];
flux = [];
for i = 1:numel(gridType)
    load(gridType{i});
    time = [time,tsave];
    flux = [flux, Wsave];
    G3D.cells.num
end


fig = figure();
hold on
figCounter = 1;
for j = [-1,0]
    subplot(1,2,figCounter);
    hold on
    figCounter = figCounter + 1;
    for i = 1:numel(gridType)
       plot(time(:,i), -flux(:,2*i +j))
    end
    
    legend(legendName)
    %axis([0, 120, 0, 3*1.1])
    xlabel('Time (s)')
    ylabel('Oil Flux ({m^3}s^{-1})')

    fig = gcf();
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 12);
end

