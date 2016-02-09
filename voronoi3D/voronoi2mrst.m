function G = voronoi2mrst(V, C, aux, name)

    assert(size(V, 2) == 3, ...
          ['Function ''%s'' is only supported in three ', ...
           'space dimensions.'], mfilename);

    % Remove auxillary cells
    C            = C(~aux)';
    cell2Node    = cumsum([1, cellfun(@numel, C)]);
    activeVertex = cell2mat(C);
    [activeVertex, ~, C] = unique(activeVertex);
    V            = V(activeVertex,:);

    % Set number of cells
    G.cells.num  = numel(cell2Node)-1;
    
    %% Find half faces
    facePos    = ones(G.cells.num+1,1);
    hf         = [];      
    hf2NodePos = [1]; 

    for i = 1:G.cells.num
        % Calculate convex hull
        hull = C(cell2Node(i)-1+convhull(V(C(cell2Node(i):cell2Node(i+1)-1),:)));
        % Merge triangle faces into polygons
        [hull, localPos] = remParFaces(V, hull);
        hf           = [hf; hull];
        hf2NodePos   = [hf2NodePos; hf2NodePos(end)-1 + localPos(2:end)];
        facePos(i+1) = numel(hf2NodePos-1);
    end
    G.cells.facePos = facePos;
    
    %% Find faces
    [nodes, nodePos, ic]  = uniqueFace(hf,hf2NodePos); 
    G.faces.nodePos = nodePos;
    G.faces.nodes   = reshape(nodes', [], 1);
    G.cells.faces   = ic;
    G.nodes.coords  = V;
    G.nodes.num     = size(V,1);  
    G.faces.num     = max(G.cells.faces);

    %% Set neighbors
    cellNo          = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
    G.faces.neighbors = zeros(G.faces.num,2);
    for i = 1:G.faces.num
        neigh = G.cells.faces==i;
        if sum(neigh)==2
            G.faces.neighbors(i,:) = cellNo(neigh);
        else
            G.faces.neighbors(i,:) = [cellNo(neigh),0];
        end
    end
    
    %% Set grid info
    G.type    = { name };
    G.griddim = 3;
end


function [faces, I] = faceSort(faces, nodePos)
    % Sort indexes in each face
    I = zeros(size(faces));
    for i = 1:size(nodePos)-1
        [faces(nodePos(i):nodePos(i+1)-1),I(nodePos(i):nodePos(i+1)-1)] ...
            = sort(faces(nodePos(i):nodePos(i+1)-1));
        I(nodePos(i):nodePos(i+1)-1) = I(nodePos(i):nodePos(i+1)-1) + nodePos(i)-1;
    end
end


function [faceNodes,nodePos, ic] = uniqueFace(hfNodes, hf2node)
    % Finds the faces given half faces.
    
    % Find all half faces with equal number of nodes
    faceSize  = diff(hf2node);
    [~,ias,~] = unique(faceSize);

    faceNodes = [];
    nodePos   = [1];
    ic = zeros(size(hf2node,1)-1,1);
    for i = 1:numel(ias) % For all half faces with equal number of nodes
        % Find the indexes of the half face nodes
        testPos  = faceSize(ias(i))==faceSize;
        fromFace = hf2node([testPos;false]);
        toFace   = hf2node([false;testPos]) - 1;
        nodeID   = arrayfun(@(l,r) (l:r), fromFace, toFace, 'un', 0)';
        nodeID   = cell2mat(nodeID);
        tempFace = reshape(hfNodes(nodeID),faceSize(ias(i)),[]);
        tempFace = tempFace';
        
        % Half faces with the same nodes are one face
        [~,ia2,ic2] = unique(sort(tempFace,2),'rows');
        temp        = tempFace(ia2,:)';
        faceNodes   = [faceNodes;temp(:)];
        ic(testPos) = ic2+numel(nodePos)-1;
        nodePos     = [nodePos; ...
                       nodePos(end)+cumsum(repmat(faceSize(ias(i)),[size(temp,2),1]))];

    end
end

function [newHull, nodePos] = remParFaces(V, hull)
    % Merge all faces in convex hull that have equal normals
    
    newHull = [];
    nodePos = [1];
    % Calculate normals
    n       = calcNormals(hull, V);
    
    while ~isempty(hull)
        % Find face normals that are equal to first face normal in stack
        parFace = n*(n(1,:)')> 1 - 50*eps;
        % Merge parallel faces
        tmp     = mergeFaces(V, hull, parFace, n(1,:)');
        nf      = numel(tmp);
        newHull = [newHull; tmp];
        nodePos = [nodePos; nodePos(end) + nf];
        % Update hull
        hull    = hull(~parFace,:);
        n       = n(~parFace,:);
    end
end


function [merged] = mergeFaces(V, H, F, n)
    % Merge faces nodes in counterclockwise direction
    
    % shif coordinate system
    merged     = unique(reshape(H(F,:),[],1));
    x0         = mean(V(merged,:));
    V          = bsxfun(@minus, V, x0);

    % Create new basis
    id         = find(F);
    basis      = [V(H(id(1),1:2),:)',n];
    basis(:,1) = basis(:,1)/norm(basis(:,1),2);
    basis(:,2) = basis(:,2)-(basis(:,1)'*basis(:,2))*basis(:,1);
    basis(:,2) = basis(:,2)/norm(basis(:,2),2);
    
    % find coordinates in new basis
    VB         = (basis\V(merged,:)')';
    
    % Sort nodes based on the angle
    VB         = VB(:,[1,2]); %Remove normal coordinate (this is zero)
    theta      = atan2(VB(:,2),VB(:,1));
    [~,i]      = sort(theta);
    merged     = merged(i);    
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

