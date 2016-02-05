function [G] = voronoi3D(pts, bound)
    % function creates 3D voronoi diagram inside bound
    
    % Calculate convex hull and remove all points outside domain.
    dt = delaunayTriangulation(bound);
    K = convexHull(dt);
    inHull = ~isnan(pointLocation(dt, pts));
    pts = pts(inHull,:);

    % Calculate boundary normals
    normals = calcNormals(K, bound);
    % Find and remove boundary planes that are the same
    remove = remParPlanes(bound(K(:,1),:), normals);
    K = K(~remove,:);
    normals = normals(~remove,:);
    
    % mirror points around planes
    for i = 1:size(K,1)
        pts = [pts; mirror(pts, bound(K(i,1),:),normals(i,:))];
    end
    
    % Generate voronoi grid
    [V, C] = voronoin(pts);

    % Find auxillary cells
    remove = false(numel(C),1);
    for i = 1:numel(C)
        avgCell = sum(V(C{i},:),1)/size(C{i},2);
        if any(isinf(avgCell)) || isnan(pointLocation(dt, avgCell));
            remove(i) = true;
        end
    end
    
    % Convert to mrst Grid-structure
    G = voronoi2mrst(V,C, remove, mfilename);

end


function r = remParPlanes(x0, n)
    r = false(size(x0,1),1);
    for i = 1:size(n,1)-1
        % test if point is on plane
        isOnPlane = (bsxfun(@minus, x0(i,:), x0(i+1:end,:))*n(i,:)') < 50*eps;
        if any(isOnPlane)
            % test normals are parallel
            isOnPlane = [false(i,1);isOnPlane];
            r(i) = any(abs(n(isOnPlane,:)*(n(i,:)'))> 1 - 50*eps);
        end
    end
end


function newPts = mirror(pts, x0, n)
    % Mirrors pts on plane
    newPts = pts - 2*bsxfun(@minus, pts, x0)*n'*n;
end