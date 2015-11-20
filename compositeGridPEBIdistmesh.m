function varargout = compositeGridPEBI(gridSizes, pdims, varargin)
        
    opt = struct('padding', 1, ...
                 'wellLines', {{}}, ...
                 'mlqtMaxLevel', 0, ...
                 'mlqtLevelSteps', -1, ...
                 'faultLines', {{}}, ...
                 'faultGridSize', -1, ...
                 'circleFactor', 0.6, ...
                 'fullFaultEdge', 1);         
    opt = merge_options(opt, varargin{:});

    assert(numel(gridSizes) == 3)

    resGridSize = gridSizes(1);
    wellGridSize = gridSizes(2);
    faultGridSize = gridSizes(3);

    faultType = [];
    wellType = [];
    priIndex = [];
    gridSpacing = [];

    
    %%
    %%Create well grid points
    if wellGridSize == -1
        wellGridSize = resGridSize/2;
    end
    
    wellPts = [];
    for i = 1:numel(opt.wellLines)
        wellLine = opt.wellLines{i};
        [newWellPts,wellSpace] = createWellGridPoints(wellLine, wellGridSize);
        np = size(newWellPts,1);
        faultType = [faultType; zeros(np,1)];
        wellType = [wellType; ones(np,1)];
        priIndex = [priIndex; (2+i)*ones(np,1)];
        gridSpacing = [gridSpacing; wellSpace];
        wellPts= [wellPts; newWellPts];
    end
    
    %%
    %%Create Fault grid points:
    faultPts = []; 
    faultCenter = [];
    faultRadius = [];
    faultToCenter = zeros(size(wellPts,1),1);
    faultPos = size(wellPts,1)+1;
    if faultGridSize ==-1
        faultGridSize = resGridSize/2;
    end
    if size(priIndex,1)<1
        maxWellPri = 2;
    else
        maxWellPri = priIndex(end);
    end

    lastFaultType = 0;
    lastCCid = 0;
    for i = 1:numel(opt.faultLines)
        fracLine = opt.faultLines{i};    
        [newFracPts, fracSpace, CC, CR,CCid] = createFracGridPoints(fracLine,...
                                                               faultGridSize,...
                                                               opt.circleFactor);

        nl = size(newFracPts,1)/2;
        if nl==0
            continue
        end
        
        newFaultType = lastFaultType+1:lastFaultType+nl;
        lastFaultType = newFaultType(end);
        faultType = [faultType; newFaultType';newFaultType'];        % 2*i-1rldecode([2*i-1; 2*i], [nl; nl])];
        wellType  = [wellType; zeros(2*nl,1)];
        priIndex = [priIndex; (maxWellPri+i)*ones(2*nl,1)];
        gridSpacing = [gridSpacing;fracSpace];
        faultPts = [faultPts;newFracPts];
        faultCenter = [faultCenter; CC];
        faultRadius = [faultRadius; CR];
        faultToCenter = [faultToCenter;lastCCid+CCid]; 
        lastCCid = faultToCenter(end)+1;
        faultPos = [faultPos; faultPos(end)+size(newFracPts,1)];
    end

    
    %% Remove fault and well conflic points
    fixedPts = [wellPts;faultPts];
    [fixedPts, removed, wellType] = removeConflictPoints(fixedPts, gridSpacing, priIndex, wellType);
    
    
    faultType = faultType(~removed);
    wellType = wellType(~removed);
    
    %%
    %% Create grid
    % set dist function
    x = pdims;
    fd = @(p) drectangle(p,0,x(1),0,x(2));
    corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];
    fixedPts = [fixedPts; corners];
    gridSpacing = [gridSpacing; zeros(size(corners,1),1)];
    priIndex = [priIndex; zeros(size(corners,1),1)];
     faultType = [faultType;zeros(size(corners,1),1)];
    wellType = [wellType;zeros(size(corners,1),1)];
    
    [Pts,t] = distmesh2d(fd, @huniform, resGridSize,[0,0;x(1),x(2)], fixedPts);
        
    G = triangleGrid(Pts, t);
    G = pebi(G);

%     %wellType = wellType + cellsContPts(G, wellPts(removed(1:size(wellPts,1)),:));
%     %label fault faces.
%     N = G.faces.neighbors + 1;
%     faultType = [0; faultType];
%     ft1 = faultType(N(:,1));
%     ft2 = faultType(N(:,2));
%     G.faces.tag = logical(ft1 == ft2 & ft1 > 0 & ft2 >0);
%     
%    %Label well cells
    %G.cells.tag = logical(wellType);

    varargout{1} = G;
    if nargout > 1
        varargout{2} = indicator;
    end
    
    figure()
     hold on
%     plot(Pts(logical(wellType),1),Pts(logical(wellType),2),'.')
%        figure()
%     hold on
%        %       hold on
%       plot(Pts(:,1),Pts(:,2),'r.')
    plot(fixedPts(:,1), fixedPts(:,2),'r.')
    plotGrid(G,'facecolor','none')
    
    n = 50;
    theta = (linspace(0,2*pi,n))';
    for i = 1:size(fixedPts,1)
        x =  fixedPts(i,1) + gridSpacing(i)*cos(theta);
        y = fixedPts(i,2) + gridSpacing(i)*sin(theta);
            plot(x,y);
        
    end
end


function [Pts, removed, wellType] = removeConflictPoints(Pts, gridSpacing, ...
                                                         priIndex, wellType)
    gridSpacing = gridSpacing*(1-1e-8); % To avoid floating point errors
    
    Ic = 1:size(Pts, 1);
    ptsToClose = Pts;
    removed = false(size(Pts, 1), 1);
    
    distance = pdist(ptsToClose)';
    dlt = distLessThan(distance, gridSpacing(Ic));
    Ic = findToClose(dlt);
    while length(Ic)>1
        %sumToClose = sumToClosePts(dlt);
        %sumToClose = sumToClose(find(sumToClose));
        %[~, Is] = sort(sumToClose,'descend');
        [~, Ii ] = sort(priIndex(Ic), 'ascend');

        removePoint = Ic(Ii(1));
        if wellType(removePoint)
            n = size(Ic,1);
            if Ii(1)>n/2;
                wellType(Ic(ceil(mod(Ii(1),(n+1)/2)))) = true;
            else
                wellType(Ic(Ii(1)+n/2)) = true;
            end
        end
        removed(removePoint) = true;    
        Ic = Ic(Ic~=removePoint);
        Ic = unique(Ic);
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

function [indexes] = findToClose(arr)
    n = length(arr);
    k = find(arr);
    [i, j] = arrToMat(k, n);
    indexes = [i;j];
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




function [Pts, removed] = enforceSufficientFaultCondition(Pts, isFault, CC, CR);
assert(size(CC,1)==size(CR,1));
assert(size(isFault,1)==size(Pts,1));
nc = size(CC,1);
np = size(Pts,1);
removed = zeros(np,1);
CRSqr = CR.^2;
for i = 1:nc
    distSqr = sum((repmat(CC(i,:),np,1)-Pts).^2,2);
    removed = removed + and(distSqr<CRSqr(i), ~isFault);
end                
removed = logical(removed);
Pts = Pts(~removed,:);
end