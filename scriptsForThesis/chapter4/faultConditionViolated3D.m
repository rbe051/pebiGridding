clear; close all
%% Generate two triangles sharing one face

pts = [0,0;0,1;1,0;0.75,0.75];
R =repmat(0.75,size(pts,1),1);
G = triangleGrid(pts);


%% Plotting

figure; hold on
plotGrid(G,'facecolor','none')

theta = linspace(0,2*pi)';
for j = 1:size(pts,1)
X = repmat(pts(j,:),100,1) + repmat(R(j),100,2).*[cos(theta), sin(theta)];
a = plot(X(:,1), X(:,2));
end