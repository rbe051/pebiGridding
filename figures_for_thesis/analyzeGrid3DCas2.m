close all; clear


fileName = {'case2fineCart', 'case2coarseCart','case2composite','case2distmesh' };
legendName = {'Fine Cartesian', 'Coarse Cartesian', 'Algorithm 1','Algorithm 2'};


figure()
hold on
for i = 1:numel(fileName)
    load(fileName{i});
    rate = zeros(numel(sol),1);
    t = zeros(numel(sol),1);
    pavg = zeros(numel(sol),1);
    for i = 1:numel(sol)
       solStep = sol(i);
       p = solStep.pressure;
       bhp = solStep.bhp;
       qS = solStep.qS;

       t(i) = solStep.time;
       rate(i) = qS;
       pavg(i) = mean(p);

    end
    plot(t(1:end)/day(), -rate(1:end)*day())
    
    %plot(t, pavg/barsa)
end
legend(legendName)
xlabel('Time (s)')
ylabel('Oil Flux ({m^3}s^{-1})')

fig = gcf();
set(findall(fig, '-property', 'FontSize'), 'FontSize', 12);

