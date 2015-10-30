function varargout = compositeGridPEBI(dims, pdims, varargin)
        
    opt = struct('padding', 1, ...
                 'wellLines', {{}}, ...
                 'wellGridSize', -1, ...
                 'fracLines', {{}}, ...
                 'fracGridSize', -1, ...
                 'circleFactor', 0.6);         
    opt = merge_options(opt, varargin{:});
             
    nx = dims(1);
    ny = dims(2);
    
    dx = pdims(1)/(nx - 1);
    dy = pdims(2)/(ny - 1);

    faultType = [];
    priIndex = [];
    gridSpacing = [];
    fracPts = [];
    %%
    %%Create Fault grid points:
    if opt.fracGridSize ==-1
        fracGridSize = sqrt(dx^2+dy^2);
    else
        fracGridSize = opt.fracGridSize;
    end
    

    lastFaultType = 0;
    for i = 1:numel(opt.fracLines)
        fracLine = opt.fracLines{i};    
        [fracPts, fracSpace] = createFracGridPoints(fracLine,...
                                                      fracGridSize,...
                                                      opt.circleFactor);

        nl = size(fracPts,1)/2;
        if nl==0
            continue
        end
        
        newFaultType = lastFaultType+1:lastFaultType+nl;
        lastFaultType = newFaultType(end);
        faultType = [faultType; newFaultType';newFaultType'];        % 2*i-1rldecode([2*i-1; 2*i], [nl; nl])];
        priIndex = [priIndex; (2+i)*ones(2*nl,1)];
        gridSpacing = [gridSpacing;fracSpace];
        fracPts = [fracPts;fracPts];
    end
    %%
    %%Create well grid points
    if opt.wellGridSize == -1
        wellGridSize = sqrt(dx^2+dy^2);
    else
        wellGridSize = opt.wellGridSize;
    end
    
    wellPts = []
    for i = 1:numel(opt.wellLines)
        wellLine = opt.wellLines{i};
        [newWellPts,wellSpace] = createWellGridPoints(wellLine, wellGridSize);
        np = size(newWellPts,1);
        faultType = [faultType; zeros(np,1)];
        priIndex = [priIndex; 2*ones(np,1)];
        gridSpacing = [gridSpacing; wellSpace];
        wellPts= [wellPts; newWellPts];
    end
    
    %%
    %% Create reservoir grid
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
    resPtsInit = [X(:), Y(:)];
    res = {};
    resGridSize = 0.5*min(dx,dy);
    for i = 1:size(resPtsInit,1)
        res = [res; mlqt(resPtsInit(i,:), wellPts, resGridSize,1, 2, dx)];
    end
    resPts = vec2mat([res{:,1}],2)
    resGridSize = [res{:,2}]';
    
    
    faultType = [faultType; zeros(size(resPts,1),1)];
    priIndex = [priIndex; ones(size(resPts, 1), 1)];
    gridSpacing = [gridSpacing; resGridSize];%(1-10^-6)*min(dx,dy)*ones(size(resPts,1),1)];
    
    
    %%
    % Put grid Points together
    Pts = [fracPts;wellPts;resPts];
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
    
    figure()
    plot(Pts(:,1),Pts(:,2),'.')
    plotGrid(G,'facecolor','none')
    varargout{1} = G;
    if nargout > 1
        varargout{2} = indicator;
    end
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
