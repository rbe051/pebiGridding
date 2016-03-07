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
        hull = C(cell2Node(i)-1+convhull(V(C(cell2Node(i):cell2Node(i+1)-1),:),'simplify',true));
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
    e = any(isnan(n),2);
    n = n(~e,:);
    hull = hull(~e,:);
%     a = patch('Vertices', V, 'faces', hull,'facecolor','y','facealpha',0.1);
%     hold on

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
    
%     a = [];
%     for i = 1:numel(nodePos)-1
%        nodes = newHull(nodePos(i):nodePos(i+1)-1);
%        a1 = patch('Vertices', V, 'faces', nodes','facecolor','y');
%        a = [a,a1];
%     end
%     hold on
%     for i = 1:numel(nodePos)-1
%        nodes = newHull(nodePos(i):nodePos(i+1)-1);
%        b = patch('Vertices', V, 'faces', nodes','facealpha',0.1);
%        for j=1:numel(nodes);
%            b2 = plot3(V(nodes(j),1),V(nodes(j),2),V(nodes(j),3),'.g','markersize',30);
%            b = [b, b2];
%        end
%        delete(b)
%     end
%     delete(a)
    
end


function [merged] = mergeFaces(V, H, F, n)
    % Merge faces nodes in counterclockwise direction
    
    % Find unique node index
    merged     = unique(reshape(H(F,:),[],1));
    % shift coordinate system
    id         = find(F);
    x0         = mean(V(H(id(1),:),:));
    VC          = bsxfun(@minus, V, x0);
    % Create new basis
    basis      = VC(H(id(1),1:2),:)';
    basis(:,1) = basis(:,1)/norm(basis(:,1),2);
    basis(:,2) = basis(:,2)-(basis(:,1)'*basis(:,2))*basis(:,1);
    basis(:,2) = basis(:,2)/norm(basis(:,2),2);
    
    % find coordinates in new basis
    VB         = (basis\VC(merged,:)')';
    
    % Sort nodes based on the angle
    theta      = atan2(VB(:,2),VB(:,1));
    [~,i]      = sort(theta);
    merged     = merged(i);
    
    % Remove colinear points
    merged = [merged;merged(1:2)];
    j = 1;
    k = 2;
    rem = false(size(merged,1),1);
    figure(2)
    a = plot(VB(:,1),VB(:,2),'.','markersize',30);
    
    while j<size(merged,1)-1
        if isColinear(V([merged(j);merged(k:k+1)],:));
           rem(k) = true;
           k = k+1;
        else
            j = k;
            k = j+1;
        end
        if  k==size(merged,1)
            break
        end

    end
    merged = merged([~rem(end-1);~rem(2:end-2);false;false],:);
    figure(1)
    b = plot3(V(merged,1),V(merged,2),V(merged,3),'.','markersize',30);
    delete(b);
    delete(a);

%     
%     
%     V = bsxfun(@plus, V, x0);
%     b = patch('Vertices', V, 'faces', merged','facealpha',0.1);
%     cent = sum(V(merged,:),1)/numel(merged);
%     b2 = plot3([cent(1),cent(1)+n(1)],[cent(2),cent(2)+n(2)], [cent(3),cent(3)+n(3)]);
%     b = [b,b2]
%     for j = 1:numel(merged)
%         b2 = plot3(V(merged(j),1), V(merged(j),2), V(merged(j),3),'.g','markersize',30);
%         b = [b,b2];
%     end
%     delete(b)
    
end
function [id] = isColinear(pts)
    id = rank(bsxfun(@minus,pts(2:end,:),pts(1,:)))<2;
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

