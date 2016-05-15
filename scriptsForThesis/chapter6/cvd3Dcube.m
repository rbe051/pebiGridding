clear; close all


%% set boundary 
boundary = [0,0,0;  ...
            1,0,0;  ...
            1,1,0;  ...
            0,1,0;  ...
            0,0,1;  ...
            1,0,1;  ...
            1,1,1;  ...
            0,1,1]; ...

%% set options
n = [100, 500, 1000];
varargin = {'storedVec', 10, 'tol' , 1e-6};

%% Generate grid

Gs   = cell(numel(n),1);
f    = cell(numel(n),1);
g    = cell(numel(n),1);
pOpt = cell(numel(n),1);
for i = 1:numel(n)
  pts = rand(n(i),3);
  [Gs{i}, pOpt{i}, f{i},g{i}] = optiVoronoi3D(pts, boundary, varargin{:});
end
%% Save grid
save('workspaceCvd3Dcube2');

%% Load grid
load('workspaceCvd3Dcube.mat')

%% Plot Grid

x = [0,1,1,0,0];
y = [0,0,1,1,0];
z = [0,0,0,0,0];

       
for i = 1:numel(Gs);
  figure(i); hold on
  G = Gs{i};
  G = computeGeometry(G);
  bndrFace = any(G.faces.neighbors==0,2); 
  outer = [G.faces.neighbors(G.faces.neighbors(:,1)==0,2); ...
           G.faces.neighbors(G.faces.neighbors(:,2)==0,1)];
  outer = unique(outer);
  oId = false(G.cells.num,1);
  oId(outer) = true;
  c = G.cells.centroids(:,1)>0.4;
  plotGrid(G, oId&c)
  plotGrid(G, ~oId&c,'facecolor','g')
  % Plot box
  plot3(x,y,z,'k')
  plot3(x,y,z+1,'k')
  for j = 1:numel(x)
    plot3([x(j),x(j)],[y(j),y(j)],[0,1],'k')
  end
  axis equal tight off
  view(-40,25)
  light('Position',[-1 -1 -1],'Style','infinite')
end

%% Save
for i = 1:numel(n)
  figure(i)
  print(strcat('../../../../master/thesis/fig/ch06/unitCubeGrid',num2str(n(i))),'-depsc')
end

%% Plot convergence
close all
for i = 1:numel(f);
  fPlt = f{i};
  gPlt = g{i};
  
  figure(i);
  semilogy(fPlt)
  axis([0,numel(fPlt),min(fPlt)*0.95,mean([max(fPlt),min(fPlt)])])
  xlabel('Number of iterations')
  ylabel('F(x)')
  
  figure(i + numel(f));
  semilogy(gPlt/gPlt(1))
  xlabel('Number of iterations')
  ylabel('||\nabla F(x_k)|| / ||\nabla F(x_0)||')
end

%% Save
for i = 1:numel(n)
  figure(i)
  print(strcat('../../../../master/thesis/fig/ch06/unitCubeGridOptFunc',num2str(n(i))),'-depsc')
  figure(i+numel(f))
  print(strcat('../../../../master/thesis/fig/ch06/unitCubeGridGradFunc',num2str(n(i))),'-depsc')
end