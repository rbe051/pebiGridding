close all
x = 0:11;
l = {[0.5,0.1; 0.5,0.9]};

figure()
hold on
Gp = compositeGridPEBI([11,11], [1, 1], 'lines', l, 'padding', 1,...
                       'fracgridsize', 0.2, 'circlefactor', 0.55);
% Add this to compositeGridPEBI to plot grid and fault nodes
% 
% plot(Pts(:,1), Pts(:,2),'bo')
% plot(left(:,1), left(:,2) ,'r.')
% plot(right(:,1), right(:,2),'r.')
                   
% Add this to compositeGrid Pebi to remove unwanted nodes
% Pts = Pts([1:53,63,64,74,75,85:end],:)
% gridSpacing = gridSpacing([1:53,63,64,74,75,85:end])
% priIndex = priIndex([1:53,63,64,74,75,85:end])
plotGrid(Gp, 'faceColor', 'none')
axis equal tight
hold on
l = l{1}
plot(l(:,1), l(:,2),'r')
axis equal off