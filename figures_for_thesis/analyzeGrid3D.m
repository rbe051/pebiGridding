close all; clear


gridType = {'fineCart', 'coarseCart','composite', 'distmesh'};
legendName = {'Fine Cartesian', 'Coarse Cartesian', ...
              'Algorithm 1', 'Algorithm 2'};
time = [];
sat = [];
flux = [];
for i = 1:numel(gridType)
    load(gridType{i});
    time = [time,tsave];
    sat = [sat, Wsave(:,2)];
    newFlux = state.wellSol.flux;
    flux = [flux,abs(sum(newFlux))];
    G3D.cells.num
end



figure()
hold on
for i = 1:numel(gridType)
   plot(time(:,i), sat(:,i)*flux(i))
end
legend(legendName)
axis([0, 120, 0, 3*1.1])
xlabel('Time (s)')
ylabel('Oil Flux ({m^3}s^{-1})')

fig = gcf();
set(findall(fig, '-property', 'FontSize'), 'FontSize', 12);

