function G = clippedPebi2D(p, bnd)

dt = delaunayTriangulation(p);
E = dt.edges();

m = size(bnd,1);



rem = true(size(E,1),1);
C = cell(size(p,1),1);
V = [];
for s = 1:size(p,1)
    NC = [E(:,2)==s, E(:,1)==s];
    bisect = find(any(NC,2));
  
    n = bsxfun(@minus, p(E(NC),:), p(s,:));
    n = bsxfun(@rdivide, n,sqrt(sum(n.^2,2)));
    x0 =bsxfun(@plus, p(E(NC),:), p(s,:))/2;
    
    symT = mat2cell(-ones(size(bnd,1),2),ones(size(bnd,1),1),2);
    
    if s == 5
      disp('s')
    end
    
    [newVertex, ~] = clipPolygon(bnd, n, x0, symT,[], bisect);
    
    if isempty(newVertex)
      rem(any(NC,2)) = false;
      continue
    end
    C{s} = [C{s}, size(V,1)+1:size(V,1)+size(newVertex,1)];
    V = [V;newVertex];
    
end

[V,IA,IC] = uniquetol(V,50*eps,'byRows',true);

G.cells.num = numel(C);
G.cells.facePos = cumsum([1; cellfun(@numel, C)]);
C2 = cellfun(@(c) [c,c(1)],C,'un', false);
faces = horzcat(C2{:})';
faces = [faces(1:end-1), faces(2:end)];
addTo = cumsum([0;ones(G.cells.num-2,1)]);

faces(G.cells.facePos(2:end-1) + addTo,:) = [];
faces = sort(IC(faces),2);

[nodes,~,faces] = unique(faces, 'rows');

G.cells.faces = faces;

G.nodes.coords = V;
G.nodes.num = size(V,1);

G.faces.num = max(faces);
G.faces.nodes = reshape(nodes',[],1);
G.faces.nodePos = (1:2:2*G.faces.num+1)';


cellNo            = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
G.faces.neighbors = zeros(G.faces.num,2);
for i = 1:G.faces.num
    neigh = G.cells.faces==i;
    if sum(neigh)==2
        G.faces.neighbors(i,[1,2]) = cellNo(neigh);
    else
        G.faces.neighbors(i,:) = [cellNo(neigh),0];
    end
end

G.griddim = 2;

G = sortEdges(G);
end


function [nodes, nodePos, faces] = uniqueFace(hf2n, hf2nPos)
  %% Find faces
  fSize   = diff(hf2nPos);     % Number of nodes of each half face
  [~,ias,~] = unique(fSize); % The unique sizes
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
end