
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
for i = 1:numel(l)
   line = l{i};
    plot(line(:, 1), line(:, 2));
end