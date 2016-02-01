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
    for i = 1:G.cells.num
       hf = [hf; convhull(V(C(nodePos(i):nodePos(i+1)-1),:))];
       facePos(i+1) = size(hf,1)+1;
    end
    G.cells.facePos = facePos;
    
    hf            = sort(hf,2);
    [hf, ia, ic]  = unique(hf, 'rows');
    G.faces.nodes = reshape(hf', [], 1);
    G.cells.faces = ic;
    G.nodes.coords = V;
        
    G.faces.num     = size(hf,1);
    %% All faces should have three edges!!
    
    
    %G.faces.nodePos = cumsum([1; repmat(2, [G.faces.num, 1])]);
    
    
    G.type = { 'pebi' };
    G.griddim = 3;

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

