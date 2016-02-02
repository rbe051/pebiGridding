function G = voronoi2mrst(V, C, aux)

    assert(size(V, 2) == 3, ...
          ['Function ''%s'' is only supported in three ', ...
           'space dimensions.'], mfilename);

    % Remove auxillary cells
    C = {C{~aux}};
    nodePos = cumsum([1, cellfun(@numel, C)]);
    activeVertex = cell2mat(C);
    [activeVertex, ~, C] = unique(activeVertex);
    V = V(activeVertex,:);

    % Set number of cells
    G.cells.num = numel(nodePos)-1;
    
    % Find half faces
    facePos = ones(G.cells.num+1,1);
    hf = [];
    %figure()
    for i = 1:G.cells.num
       hf = [hf; C(nodePos(i)-1+convhull(V(C(nodePos(i):nodePos(i+1)-1),:)))];
       facePos(i+1) = size(hf,1)+1;
       %patch('Vertices',V,'Faces',hf(facePos(i):facePos(i+1)-1,:),'FaceColor',[0,0,i/numel(C)],'FaceAlpha',0.5) 
    end
    G.cells.facePos = facePos;
    
    [hf, I]      = sort(hf,2);
    [hf, ~, ic]  = unique(hf, 'rows');
    G.faces.nodes = reshape(hf', [], 1);
    G.cells.faces = ic;
    G.nodes.coords = V;
    G.nodes.num    = size(V,1);  
    G.faces.num     = size(hf,1);

    G.faces.nodePos   = cumsum([1;repmat(3, [G.faces.num, 1])]);
    cellNo            = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
    
    G.faces.neighbors = zeros(G.faces.num,2);
    hold on
    for i = 1:G.faces.num
        neigh = G.cells.faces==i;
        find(neigh)
        G.nodes.coords(G.faces.nodes(G.faces.nodePos(i):G.faces.nodePos(i+1)-1),:)
        h1 = plot3(G.nodes.coords(G.faces.nodes(G.faces.nodePos(i):G.faces.nodePos(i+1)-1),1), ...
              G.nodes.coords(G.faces.nodes(G.faces.nodePos(i):G.faces.nodePos(i+1)-1),2), ...
              G.nodes.coords(G.faces.nodes(G.faces.nodePos(i):G.faces.nodePos(i+1)-1),3), ...
              'g.', 'markersize', 30);
        if sum(neigh)==2
            G.faces.neighbors(i,:) = cellNo(neigh);
        else
            G.faces.neighbors(i,:) = [cellNo(neigh),0];
        end
        G.faces.neighbors
        delete(h1)
    end
    
    
    G.type = { 'pebi' };
    G.griddim = 3;
    figure()
    plotGrid(G,'FaceColor',[0,0,1/20],'FaceAlpha',0.5)

end


function [ID] = remEqEdges(ID, n)
% Removes all edges that whose neighboor faces has equal normals

    % find normals that are equal
    for i = size(n,1)-1
        areEqual = all(abs(bsxfun(@minus, n(i,:), n(i+1:end,:))) < 50*eps);
        if any(areEqual)
           rem = all(ismember(ID(areEqual,:),ID(k,:)),1); 
        end
    end


end

