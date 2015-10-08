function varargout = compositeGridPEBI(dims, pdims, varargin)
    nx = dims(1);
    ny = dims(2);
    
    dx = pdims(1)/(nx - 1);
    dy = pdims(2)/(ny - 1);

    
    opt = struct('padding', 1, ...
                 'lines', {{}}, ...
                 'circleFactor', 0.6);
         
    opt = merge_options(opt, varargin{:});

    circleFactor = opt.circleFactor;
    assert(0.5<circleFactor && circleFactor < 1)
    
    vx = 0:dx:pdims(1);
    vy = 0:dy:pdims(2);
    [X, Y] = meshgrid(vx, vy);


    [ii, jj] = meshgrid(1:nx, 1:ny);

    nedge = opt.padding;
    exterior = (ii <= nedge | ii > nx - nedge) | ...
               (jj <= nedge | jj > ny - nedge);

    interior = ~exterior;

    X(interior) = X(interior);
    Y(interior) = Y(interior);

    Pts = [X(:), Y(:)];
    Pts0 = Pts;

    Pts = [X(:), Y(:)];
    Pts0 = Pts;
    
    priIndex = ones(size(Pts, 1), 1);
    
    faultType = zeros(size(Pts, 1), 1);
    for i = 1:numel(opt.lines)
        l = opt.lines{i};
        
        assert(all(size(l) == [2 2]));
        p1 = l(1, :);
        p2 = l(2, :);
        
        v = p2 - p1;
        
        dists = norm(v, 2)/norm([dx, dy], 2);
        dists = 2.0*max(ceil(dists), 2);
        
        % Expand into a vector
        l = p1;
        for j = 1:dists
            l = [l; p1 + v*j/dists];
        end
        
        % Create left vector and right vector(to contain fracture points).
        left = zeros(dists-1, 2);
        right = zeros(dists-1, 2);

        for j=1:dists-1
            line_length = norm(l(j+1,:)-l(j,:), 2);   %||p_(j+1) - p_j||
            n1 = (l(j+1,:)-l(j,:))/line_length;       %Unit vector
            n2 = [-n1(2), n1(1)];                     %Unit normal
            fracture_distance = 2*line_length/2*sqrt(4*circleFactor^2 -1)*n2;  %% OBS! Should make sure fracture points are equaly spaced.
            left(j,:) = l(j,:) + line_length/2*n1 + 0.5*fracture_distance;
            right(j,:) = l(j,:) + line_length/2*n1 - 0.5*fracture_distance;
        end

        nl = size(left, 1);
        [Pts, ~, removed] = replacePointsByHull(Pts, [left; right]);
        faultType = faultType(~removed);
        priIndex = priIndex(~removed);
        faultType = [faultType; rldecode([1; 2], [nl; nl])];
        priIndex = [priIndex; 10*ones(nl,1)];
    end
    
    
    Tri = delaunayTriangulation(Pts);
    
    G = triangleGrid(Pts, Tri.ConnectivityList);

    G = pebi(G);

    N = G.faces.neighbors + 1;
    faultType = [0; faultType];

    ft1 = faultType(N(:, 1));
    ft2 = faultType(N(:, 2));
    G.faces.tag = double(ft1 ~= ft2 & ft1 > 0 & ft2 > 0);
    indicator = ~ismember(Pts, Pts0, 'rows');

    varargout{1} = G;
    if nargout > 1
        varargout{2} = indicator;
    end
end


function [Pts, removed] = removeConflictPoints(Pts, gridSpacing, priIndex)
    distances = squareform(pdist(Pts));
    rm = sparse(distLessThan(distance, gridSpacing));
    One = ones(size(rm,1), 1);
    sumShort = rm*One;
    
end

function [arr] = distLessThan(distance, b)
    m = length(distance);
    k = 1:m;
    j = ceil((2*m-1)/2 - 0.5*sqrt((2*m-1)^2 - 8*k));
    i = k + j - (j-1).*(m-j/2);
    arr = distance < min(b(i), b(j));
end

function [dists] = getDistance(distance, i, j)
    %Returns the distance between point i and j
    m = 0.5*(sqrt(8*length(distance) + 1) + 1);
    assert((0<i && i<=m) && (0<j &&j<=m))
    
    if j== i
        dists = 0;
    elseif i < j
        dists = distance((i-1)*(m - 0.5*i) + j - i);
    else
        dists = distance((j-1)*(m - 0.5*j) + i - j);
    end
end


function [Pts, isNew, removed] = replacePointsByHull(Pts, P_target)
    Tri = delaunayTriangulation(P_target);
    keep = isnan(Tri.pointLocation(Pts));
    
    Pts = [Pts(keep, :); P_target];
    isNew = true(size(Pts, 1), 1);
    isNew(1:sum(keep)) = false;
    
    removed = ~keep;
end
