
%% Fault honoring pebi / triangle grids

x = 0.2:0.05:0.8;
l = {[0.65,0.1;0.65,0.926],...
     [0.2,0.175; 0.875,0.875], ...
     [0.2,0.925; 0.9,0.125], ...
     [0.45,0.15; 0.83, 0.35]};

figure()
hold on
Gp = compositeGridPEBI([30,30], [1, 1], 'lines', l, 'padding', 1);

plotGrid(Gp, 'faceColor', 'none')
axis equal tight
hold on
for i = 1:numel(l)
   line = l{i};
    %plot(line(:, 1), line(:, 2));
end