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
