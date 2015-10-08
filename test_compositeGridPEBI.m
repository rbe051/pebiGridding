
%% Fault honoring pebi / triangle grids
%close all
l = [0.50, 0.50; 0.8, 1.1];

Gp = compositeGridPEBI([15, 30], [1, 2], 'lines', {l});

figure;
plotGrid(Gp)
axis equal tight off
hold on
plot(l(:, 1), l(:, 2));
