close all
x = 0:11;
l = {[0.5,0.1; 0.5,0.9]};

figure()
hold on
Gp = compositeGridPEBI([11,11], [1, 1], 'lines', l, 'padding', 1,...
                       'fracgridsize', 0.2, 'circlefactor', 0.55);

plotGrid(Gp, 'faceColor', 'none')
axis equal tight
hold on
l = l{1}
plot(l(:,1), l(:,2),'r')
axis equal off