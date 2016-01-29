function [G] = voronoi3D(pts, bound)
    % function creates 3D voronoi diagram inside bound
    
    % Calculate convex hull and remove all points outside domain.
    dt = delaunayTriangulation(bound);
    K = convexHull(dt);
    inHull = ~isnan(pointLocation(dt, pts));
    pts = pts(inHull,:);

    % Calculate boundary normals
    normals = calcNormals(K, bound);
    % Find boundary planes that equal
    remove = remEqPlanes(bound(K(:,1),:), normals);
    K = K(~remove,:);
    normals = normals(~remove,:);
    
    % mirror points around planes
    newPts = [];
    for i = 1:size(K,1)
        newPts = [newPts; mirror(pts, bound(K(i,1),:),normals(i,:))];
    end
    
    figure
    hold on
    plot3(bound(:,1), bound(:,2), bound(:,3),'.')
    plot3(pts(:,1),pts(:,2),pts(:,3),'m.')
    pts = [pts;newPts];
    plot3(pts(:,1),pts(:,2),pts(:,3),'ro')
    
    [V, C] = voronoin(pts);
    

    for i = 1:numel(C)
       VertCell = V(C{i},:); 
       if all(C{i}~=1)
            KVert = convhulln(VertCell);
            patch('Vertices',VertCell,'Faces',KVert,'FaceColor',[0,0,i/numel(C)],'FaceAlpha',0.5) 
       end
    end
    G = [];
end

function n = calcNormals(ID, pts)
    % Calculates the normals of the triangles pts(ID)
    
    tri = pts(ID',:);
    x0  = tri(1:3:end,:);
    x1  = tri(2:3:end,:);
    x2  = tri(3:3:end,:);
    n   = cross(x1 - x0, x2 - x0);
    n   = bsxfun(@rdivide, n, sqrt(sum(n.^2,2)));
end

function r = remEqPlanes(x0, n)
    r = false(size(x0,1),1);
    for i = 1:size(n,1)-1
        % test if point is on plane
        isOnPlane = abs(bsxfun(@minus, x0(i,:), x0(i+1:end,:))*n(i,:)') < 50*eps;
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