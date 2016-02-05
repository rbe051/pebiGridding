
clc
close all


Pts = [0,0; 1,0; 0,1; -2,0; 0,-2; 1,1];
gridSpacing = 1.01*[1, 1, 1, 1, 1, 1]';
priIndex = [-1, 1, 1, 1, 1, 2]';

[pts, rem] = removeConflictPoints(Pts, gridSpacing, priIndex);

figure()
hold on
plot(Pts(:,1), Pts(:,2), 'o')
plot(pts(:,1), pts(:,2), '.')

