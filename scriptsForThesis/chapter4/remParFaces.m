function [newHull, nodePos] = remParFaces(V, hull)
    % Merge all faces in convex hull that have equal normals
    
    newHull = [];
    nodePos = [1];
    % Calculate normals
    n       = calcNormals(hull, V);

% This might be necesary
%    e = any(isnan(n),2);
%    n = n(~e,:);
%    hull = hull(~e,:);

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
    VC         = bsxfun(@minus, V, x0);
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