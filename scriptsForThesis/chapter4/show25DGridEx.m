clear;close all

%% Set fault
f = {[0.0,0.5;1,0.5]};

%% create fault
fGs = 0.1;
F = createFaultGridPoints(f,fGs);

%% Create reservoir grid points
% set domain function
x = [1,1];
rectangle = [0,0; x(1),x(1)];   
fd = @(p) drectangle(p, 0, x(1), 0, x(2));
% Set fixed points
corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];

% Add wells an faults as fixed points
fixedPts = [F.f.pts; corners];

[Pts,~,sorting] = distmesh2d(fd, @huniform, fGs, rectangle, fixedPts);

%% Create Grid
bound = [0,0;0,x(2);x(1),x(2);x(1),0];
G = clippedPebi2D(Pts,bound);

%% Create 2.5D grid
n = 4;
h = 0.5;
slope = 1.5;
Gl = makeLayeredGrid(G,n);
Gl.nodes.coords(:,3) = Gl.nodes.coords(:,3)*h/4;
Gd = Gl;
Gd.nodes.coords(:,2) = Gl.nodes.coords(:,2)+...
  slope*Gl.nodes.coords(:,end).*Gl.nodes.coords(:,2).*(x(2) - Gl.nodes.coords(:,2));

lS = [G.nodes.coords,zeros(G.nodes.num,1)];
lE = [lS(:,1),lS(:,2)+ slope*h*G.nodes.coords(:,2).*(x(2) - G.nodes.coords(:,2)), ...
      h*ones(G.nodes.num,1)];

%% Plot
figure(); hold on
plotGrid(G,'facecolor','none')
axis equal off tight
for i = 1:size(lines,1)
  plot3([lS(i,1),lE(i,1)],[lS(i,2),lE(i,2)], [lS(i,3),lE(i,3)],'k')
end

figure()
plotGrid(Gd)
axis equal off tight
view(70,30)
light('position',[1.2,0.5,-1],'style','local')
light('position',[1.5,-0.5,-0.1],'style','local')

%% Save
print('../../../../master/thesis/fig/ch04/show25DgridEx','-depsc');