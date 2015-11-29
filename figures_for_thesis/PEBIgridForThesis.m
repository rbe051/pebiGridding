clc; close all


%% Show error in pebi function
x = 1:3;
[X,Y] = meshgrid(x,x);
pts = [X(:), Y(:)];
pts(2,:) = [1.5,2];
pts(5,:) = [1.5,1.48];
%pts(6,:) = [1.2,1.9];
%pts(10,:) = [1.21, 1.99];
%pts = randn(10,2);
Gt = triangleGrid(pts);

Gp = pebi(Gt);

plotGrid(Gp, 'facecolor','none');
hold on
plot(pts(:,1), pts(:,2), 'o')
figure
plotGrid(Gt,'facecolor','none');
hold on
plot(pts(:,1), pts(:,2), 'o')