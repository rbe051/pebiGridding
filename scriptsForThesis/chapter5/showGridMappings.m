clear; close all

rng(3);
pts = rand(6,2);

G = triangleGrid(pts);
G = computeGeometry(G);
plotGrid(G,'facecolor',[0.9,0.9,0.9])
hold on


cells = num2str(1:G.cells.num);
cells(ismember(cells,' ,.:;!')) = [];
faces = num2str(1:G.faces.num);
faces(ismember(faces,' ,.:;!')) = [];
G.cells.facePos
G.cells.faces

for i = 1:G.cells.num
  text(G.cells.centroids(i,1),G.cells.centroids(i,2),num2str(i))
end
for i = 1:G.faces.num
  text(G.faces.centroids(i,1),G.faces.centroids(i,2),num2str(i))
end


axis equal off tight

%% Save
set(0,'defaulttextinterpreter','none')
%print('../../../../master/thesis/fig/ch05/inkscape/gridMappings','-dsvg')