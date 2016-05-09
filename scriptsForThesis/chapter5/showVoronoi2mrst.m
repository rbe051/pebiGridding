clear; close all
%addpath ../../voronoi3D/

%% Create voronoi Grid:
n   = 5;
[X,Y,Z] = ndgrid(linspace(0,1,n));
pts = [X(:),Y(:),Z(:)];
pts(1:2:end,1) = pts(1:2:end,1) + 0.1;
[V,C] = voronoin(pts);

% Remove infinity cells
rem = cellfun(@(c) any(isinf(V(c,1))), C);
C   = C(~rem);
V   = V(2:end,:);
C   = cellfun(@(c) c-1, C,'un',false);


%% Plot Grid
for i = 1:numel(C)
  hull = convhull(V(C{i},:));
  patch('Vertices', V(C{i},:), 'faces',hull,'facecolor','y')
end
set(gca,'zdir','reverse')
view(140,30)
axis equal tight off
light('position',[0,1,0])


%% Create the convex hull of a cell
cNum = numel(C);        % Number of cells
hf2n    = [];           % Map from half-face nodes to nodes
hf2nPos = 1;            % Map from half-faces to half-face nodes
facePos = ones(cNum,1);

for i = 1:cNum;
  hull         = convhull(V(C{i},:));
  [hull, localPos] = remParFaces(V(C{i},:), hull);
  hf2n         = [hf2n; C{i}(hull)'];
  hf2nPos      = [hf2nPos; (hf2nPos(end) - 1 + localPos(2:end))];
  facePos(i+1) = numel(hf2nPos-1);
end

%% Find faces
fSize   = diff(hf2nPos);     % Number of nodes of each half face
[~,ias,ics] = unique(fSize); % The unique sizes
nodes   = [];
nodePos = 1;
faces   = zeros(size(hf2nPos,1)-1,1);

for i = 1:numel(ias) 
  % Find the indexes of the half-face nodes 
  testPos  = fSize(ias(i))==fSize;
  from     = hf2nPos([testPos;false]);
  to       = hf2nPos([false;testPos]) - 1;
  map      = mcolon(from, to);
  nTmp     = reshape(hf2n(map),fSize(ias(i)),[])';
  % Half faces with the same nodes are one face
  [~,ia,ic]= unique(sort(nTmp,2), 'rows');
  newNodes = nTmp(ia,:)';
  nodes    = [nodes; newNodes(:)];
  faces(testPos) = ic+numel(nodePos)-1;
  locPos   = cumsum(repmat(fSize(ias(i)),[size(newNodes,2),1]));
  nodePos  = [nodePos; ...
              nodePos(end) + locPos];
end
fNum = numel(nodePos)-1;

%% Find neighbors

cellNo    = rldecode(1:cNum, diff(facePos), 2).';
neighbors = zeros(fNum,2);
for i = 1:fNum
  neigh = faces==i;
  if sum(neigh)==2
     neighbors(i,[1,2]) = cellNo(neigh);
  else
     neighbors(i,:) = [cellNo(neigh),0];
  end
end

%% create grid
G.cells.num       = cNum;
G.cells.facePos   = facePos;
G.cells.faces     = faces;

G.faces.nodePos   = nodePos;
G.faces.num       = fNum;
G.faces.nodes     = nodes;

G.nodes.num       = size(V,1);
G.nodes.coords    = V;
G.faces.neighbors = neighbors;

G.griddim         = 3;

%% Plot grid
figure()
plotGrid(G)
view(140,30)
axis equal tight off
light('position',[0,1,0])

%% Save to file
% figure(1)
% print('../../../../master/thesis/fig/ch04/voronoi2mrstPatch','-depsc')
% figure(2)
%print('../../../../master/thesis/fig/ch04/voronoi2mrstGrid','-dpng')
 

    
    
