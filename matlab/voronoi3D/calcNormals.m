function n = calcNormals(ID, pts)
    % Calculates the normals of the triangles pts(ID)
    % Arguments:
    %   ID      A triangulation index, e.g. ID = convhull(pts)
    %   pts     A nx3 array of the coordinats of the vertexes
    % returns:
    %   n       Normals of each triangle
    tri = pts(ID',:);
    x0  = tri(1:3:end,:);
    x1  = tri(2:3:end,:);
    x2  = tri(3:3:end,:);
    n   = cross(x1 - x0, x2 - x0);
    n   = bsxfun(@rdivide, n, sqrt(sum(n.^2,2)));
end