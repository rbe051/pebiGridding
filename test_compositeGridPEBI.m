
%% Fault honoring pebi / triangle grids
close all
l = {[0.50, 0.50; 0.8, 1.1], ...
      [0.5, 0.2; 0.5, 1.5]};

figure()
hold on
Gp = compositeGridPEBI([15, 30], [1, 2], 'lines', l);

plotGrid(Gp, 'faceColor', 'none')
axis equal tight
hold on
plot(l(:, 1), l(:, 2));
