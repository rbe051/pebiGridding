function G = voronoi2mrst(V, C, aux)

    assert(size(V, 2) == 3, ...
          ['Function ''%s'' is only supported in three ', ...
           'space dimensions.'], mfilename);

    % Remove auxillary cells
    C = {C{~aux}};
    cell2Node = cumsum([1, cellfun(@numel, C)]);
    activeVertex = cell2mat(C);
    [activeVertex, ~, C] = unique(activeVertex);
    V = V(activeVertex,:);

    % Set number of cells
    G.cells.num = numel(cell2Node)-1;
    
    % Find half faces
    facePos = ones(G.cells.num+1,1);
    hf = [];      
    hf2Node = [1]; 
    %figure()
    for i = 1:G.cells.num
        hull = C(cell2Node(i)-1+convhull(V(C(cell2Node(i):cell2Node(i+1)-1),:))); 
        [hull, localPos] = remParFaces(V, hull);
        hf = [hf; hull];
        hf2Node = [hf2Node; hf2Node(end)-1 + localPos(2:end)];
        facePos(i+1) = numel(hf2Node);
        for j = facePos(i):facePos(i+1)-1
        b = plot3(V(hf(hf2Node(j):hf2Node(j+1)-1),1), ...
                  V(hf(hf2Node(j):hf2Node(j+1)-1),2), ...
                  V(hf(hf2Node(j):hf2Node(j+1)-1),3), ...
                  'mo','markersize',30);
            delete(b);
        end
       %patch('Vertices',V,'Faces',hf(facePos(i):facePos(i+1)-1,:),'FaceColor',[0,0,i/numel(C)],'FaceAlpha',0.5) 
    end
    G.cells.facePos = facePos;
    
  
    %[hf, I]      = faceSort(hf,hf2Node); %I believe I sort in line 10, so this is unnessesary
    [hf, nodePos, ic]  = uniqueFace(hf,hf2Node,V);
    
    G.faces.nodePos = nodePos;
    for j = 1:numel(facePos)-1
        faces = ic(facePos(j):facePos(j+1)-1)
        for i = 1:numel(faces)
            b = plot3(V(hf(nodePos(faces(i)):nodePos(faces(i)+1)-1),1), ...
                      V(hf(nodePos(faces(i)):nodePos(faces(i)+1)-1),2), ...
                      V(hf(nodePos(faces(i)):nodePos(faces(i)+1)-1),3), ...
                        'mo','markersize',30);
            delete(b);
        end
    end
    
    G.faces.nodes = reshape(hf', [], 1);
    G.cells.faces = ic;
    G.nodes.coords = V;
    G.nodes.num   = size(V,1);  
    G.faces.num   = max(G.cells.faces);

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


function [faces, I] = faceSort(faces, nodePos)
    I = zeros(size(faces));
    for i = 1:size(nodePos)-1
        [faces(nodePos(i):nodePos(i+1)-1),I(nodePos(i):nodePos(i+1)-1)] ...
            = sort(faces(nodePos(i):nodePos(i+1)-1));
        I(nodePos(i):nodePos(i+1)-1) = I(nodePos(i):nodePos(i+1)-1) + nodePos(i)-1;
    end
end


function [newNodes,nodePos, ic] = uniqueFace(nodes, hf2node, V)
    faceSize = diff(hf2node);
    [~,ias,~] = unique(faceSize);

    newNodes = [];
    nodePos  = [1];
    ic = zeros(size(hf2node,1)-1,1);
    for i = 1:numel(ias)
        testPos = faceSize(ias(i))==faceSize;
        fromFace = hf2node([testPos;false]);
        toFace   = hf2node([false;testPos]) - 1;
        faceID  = arrayfun(@(l,r) (l:r), fromFace, toFace, 'un', 0)';
        faceID  = cell2mat(faceID);

        testFace = reshape(nodes(faceID),faceSize(ias(i)),[]);
        
        testFace = testFace';
        [aaa,ia2,ic2] = unique(sort(testFace,2),'rows');
        temp = testFace(ia2,:)';
        newNodes = [newNodes;temp(:)];
        ic(testPos) = ic2+numel(nodePos)-1;
        nodePos = [nodePos; ...
                   nodePos(end)+cumsum(repmat(faceSize(ias(i)),[size(temp,2),1]))];

    end
    
    %% plot for testing    
    for j = 1:numel(nodePos)-1
        b = plot3(V(newNodes(nodePos(j):nodePos(j+1)-1),1), ...
                  V(newNodes(nodePos(j):nodePos(j+1)-1),2), ...
                  V(newNodes(nodePos(j):nodePos(j+1)-1),3), ...
                  'mo','markersize',30);
            delete(b);
        end
end

function [newHull, nodePos] = remParFaces(V, hull)
    newHull = [];
    nodePos = [1];
    % Calculate normals
    n = calcNormals(hull, V);
    while ~isempty(hull)
        % Find face normals that are equal to top face normal
        parFace = n*(n(1,:)')> 1 - 50*eps;
        % Merge parallel faces
        tmp = mergeFaces(V, hull, parFace, n(1,:)');

        nf  = numel(tmp);
        hold on
        for i = 1:numel(tmp)
            b= plot3(V(tmp(i),1), V(tmp(i),2), V(tmp(i),3),'g.','markersize', 30);
            delete(b)
        end
        newHull = [newHull; tmp];
        nodePos = [nodePos; nodePos(end) + nf]
        % Update hull
        hull = hull(~parFace,:);
        n    = n(~parFace,:);
    end
end


function [merged] = mergeFaces(V, H, F, n)
    % Merge faces in counterclockwise direction
    
    % shif coordinate system
    merged = unique(reshape(H(F,:),[],1));
    x0 = mean(V(merged,:));
    V = bsxfun(@minus, V, x0);

    % Create new basis
    id = find(F);
    basis = [V(H(id(1),1:2),:)',n];
    basis(:,1) = basis(:,1)/norm(basis(:,1),2);
    basis(:,2) = basis(:,2)-(basis(:,1)'*basis(:,2))*basis(:,1);
    basis(:,2) = basis(:,2)/norm(basis(:,2),2);
    
    % find new coordinates
    VB = (basis\V(merged,:)')';
    %Remove normal coordinate (this is zero)
    VB = VB(:,[1,2]);
    theta = atan2(VB(:,2),VB(:,1));
    [~,i] = sort(theta);
    merged = merged(i);    
end

function [ID] = remEqEdges(V, n)
% Removes all edges that whose neighboor faces has equal normals

    % find normals that are equal
    for i = size(n,1)-1
        areEqual = all(abs(bsxfun(@minus, n(i,:), n(i+1:end,:))) < 50*eps);
        if any(areEqual)
           rem = all(ismember(ID(areEqual,:),ID(k,:)),1); 
        end
    end


end

