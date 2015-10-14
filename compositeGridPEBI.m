function varargout = compositeGridPEBI(dims, pdims, varargin)
    nx = dims(1);
    ny = dims(2);
    
    dx = pdims(1)/(nx - 1);
    dy = pdims(2)/(ny - 1);

    
    opt = struct('padding', 1, ...
                 'lines', {{}}, ...
                 'fracGridSize', -1, ...
                 'circleFactor', 0.6);
         
    opt = merge_options(opt, varargin{:});

    circleFactor = opt.circleFactor;
    fracDs = opt.fracGridSize;
    if fracDs == -1
        fracDs = min(dx,dy);
    end
    assert(0.5<circleFactor && circleFactor < 1)
    assert(0<fracDs)
    
    Pts = [];
    priIndex = [];
    gridSpacing = [];
    faultType = [];
    lastFaultType = 0;
    for i = 1:numel(opt.lines)
        fracLine = opt.lines{i};
        
        assert(all(size(fracLine,2) == 2));
        assert(all(size(fracLine,1) > 1));
        
        
        [fracLine, fracDs] = eqInterpret(fracLine, fracDs);
        %This is equidistante if you follow the line described by fracLine,
        %but the new points may be a bit to close if the fracture has
        %sharp corners and/or you upsample the line.
        
        numOfFracLine = size(fracLine,1);
        
        % Create left vector and right vector(to contain fracture points).
        left = zeros(numOfFracLine-1, 2);
        right = zeros(numOfFracLine-1, 2);
        ptGrSp = zeros(numOfFracLine-1, 1);
        for j=1:numOfFracLine-1
            lineLength = norm(fracLine(j+1,:) - fracLine(j,:));
            %Because the line lengths are not completely uniform we have to
            %calculate this for each segment.
            fractureRadius = lineLength/2*sqrt(4*circleFactor^2 -1);
            n1 = (fracLine(j+1,:)-fracLine(j,:))/lineLength;   %Unit vector
            n2 = [-n1(2), n1(1)];                              %Unit normal
            left(j,:) = fracLine(j,:) + lineLength/2*n1 + fractureRadius*n2;
            right(j,:) = fracLine(j,:) + lineLength/2*n1 - fractureRadius*n2;
            ptGrSp(j) = (2-10^-6*fracDs)*fractureRadius;
        end

        nl = (numOfFracLine-1);
        Pts = [Pts;left;right];
        newFaultType = lastFaultType+1:lastFaultType+nl;
        faultType = [faultType; newFaultType';newFaultType'];        % 2*i-1rldecode([2*i-1; 2*i], [nl; nl])];
        lastFaultType = faultType(end);
        priIndex = [priIndex; (2+i)*ones(2*nl,1)];
        gridSpacing = [gridSpacing; ptGrSp; ptGrSp];
    end
    
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

    resPts = [X(:), Y(:)];
    Pts = [Pts;resPts];
    faultType = [faultType; zeros(size(resPts,1),1)];
    priIndex = [priIndex; ones(size(Pts, 1), 1)];
    gridSpacing = [gridSpacing; (1-10^-6)*min(dx,dy)*ones(size(Pts,1),1)];
    

    [Pts, removed] = removeConflictPoints(Pts, gridSpacing, priIndex);
    faultType = faultType(~removed);
    Tri = delaunayTriangulation(Pts);
    G = triangleGrid(Pts, Tri.ConnectivityList);
    G = pebi(G);
    
    N = G.faces.neighbors + 1;
    faultType = [0; faultType];
    ft1 = faultType(N(:,1));
    ft2 = faultType(N(:,2));
    G.faces.tag = double(ft1 == ft2 & ft1 > 0 & ft2 >0);

    varargout{1} = G;
    if nargout > 1
        varargout{2} = indicator;
    end
end

function [newPoints, dt] = eqInterpret(path, dt)
    linesDist = sqrt(sum(diff(path,[],1).^2,2));
    linesDist = [0; linesDist]; % add the starting point
    cumDist = cumsum(linesDist);
    dt = cumDist(end)/ceil(cumDist(end)/dt);
    newPointsLoc = 0:dt:cumDist(end);
    if length(newPointsLoc)<2 % Fracture to short
        error('Fracture grid size larger than fracture')
    end
        
    newPoints = interp1(cumDist, path, newPointsLoc);
      
end


function [Pts, removed] = removeConflictPoints(Pts, gridSpacing, priIndex)
    Ic = 1:size(Pts, 1);
    ptsToClose = Pts;
    removed = zeros(size(Pts, 1), 1);
    
    distance = pdist(ptsToClose)';
    dlt = distLessThan(distance, gridSpacing(Ic));
    Ic = findToClose(dlt);
    
    while length(Ic)>1
        %sumToClose = sumToClosePts(dlt);
        %sumToClose = sumToClose(find(sumToClose));
        %[~, Is] = sort(sumToClose,'descend');
        [~, Ii ] = sort(priIndex(Ic), 'ascend');

        removePoint = Ic(Ii(1));
        removed(removePoint) = 1;        
        Ic = Ic(Ic~=removePoint);
        ptsToClose = Pts(Ic,:);
        
        if size(ptsToClose,1) ==1
            continue
        end
        distance = pdist(ptsToClose)';
        dlt = distLessThan(distance, gridSpacing(Ic));
        Ic = Ic(findToClose(dlt));
    end
    Pts = Pts(~removed,:);
end


function [arr] = distLessThan(distance, b)
    n = length(distance);
    [i,j] = arrToMat(1:n, n);
    arr = distance < max(b(i), b(j));
end

function [pts] = sumToClosePts(arr)
    n = length(arr);
    m = ceil(sqrt(2*n)); % = 0.5 + 0.5sqrt(1+8n) for n > 0
    pts = zeros(m,1);
    for i = 1:m
        k1 = matToArr(i+1:m, i, m);
        k2 = matToArr(i,1:i-1, m);
        pts(i) = sum(arr(k1)) + sum(arr(k2)); 
    end
end

function indexes = findToClose(arr)
    n = length(arr);
    k = find(arr);
    [i, j] = arrToMat(k, n);
    indexes = unique([i ; j]);
end

function [k] = matToArr(i,j, m)
    assert(all(abs(j)) && all(abs(i-j)) && all(abs(1+m-i)));
    k = 1 + (j-1)*m - (j-1).*j/2 + i-j - 1;
end

function [i, j] = arrToMat(k, n)
    m = ceil(sqrt(2*n));
    j = ceil((2*m-1)/2 - 0.5*sqrt((2*m-1)^2 - 8*k));
    i = k + j - (j-1).*(m-j/2);
end


function [Pts, isNew, removed] = replacePointsByHull(Pts, P_target)
    Tri = delaunayTriangulation(P_target);
    keep = isnan(Tri.pointLocation(Pts));
    
    Pts = [Pts(keep, :); P_target];
    isNew = true(size(Pts, 1), 1);
    isNew(1:sum(keep)) = false;
    
    removed = ~keep;
end
