function varargout = compositeGridPEBI(resGridSize, pdims, varargin)
        
    opt = struct('wellLines', {{}}, ...
                 'wellGridFactor', -1, ...
                 'mlqtMaxLevel', 0, ...
                 'mlqtLevelSteps', -1, ...
                 'faultLines', {{}}, ...
                 'faultGridFactor', -1, ...
                 'circleFactor', 0.6, ...
                 'fullFaultEdge', 1,...
                 'priOrder', []);         
    
    %% Set options
    opt = merge_options(opt, varargin{:});
    % Set grid sizes
    wellGridFactor = opt.wellGridFactor;
    if wellGridFactor < 0
        wellGridFactor = 0.5;
    end
    faultGridFactor = opt.faultGridFactor;
    if faultGridFactor < 0
        faultGridFactor = 0.5;
    end
    wellGridSize = resGridSize*wellGridFactor;
    faultGridSize = resGridSize*faultGridFactor;
    mlqtMaxLevel = opt.mlqtMaxLevel;
    mlqtLevelSteps = opt.mlqtLevelSteps;
    
    % Load faults and Wells
    faultLines = opt.faultLines;
    wellLines = opt.wellLines;
    nFault = numel(faultLines);
    nWell = numel(wellLines);
    linesToGrid = [wellLines, faultLines{:}];
    isWell = [true(nWell,1); false(nFault,1)];
    
    % Set priority index
    priOrder = opt.priOrder;
    if isempty(priOrder)
        priOrder = 1:nFault+nWell;
    else
        assert(numel(priOrder) == nFault + nWell);
    end
    % Sort faults and wells in priority
    linesToGrid = linesToGrid(priOrder);
    isWell = isWell(priOrder);
    


    %% Initialize variables.
    faultType = [];
    wellType = logical([]);
    priIndex = [];
    gridSpacing = [];

    Pts = [];
    faultCenter = [];
    faultRadius = [];
    faultToCenter = [];
    faultPos = [];


    %% Place fault and well points
    lastFaultType = 0;
    lastCCid = 0;
    for i = 1:nFault + nWell
        if isWell(i)
            wellLine = linesToGrid{i};
            [wellPts,wellSpace] = createWellGridPoints(wellLine, wellGridSize);
            np = size(wellPts,1);
            faultType = [faultType; zeros(np,1)];
            wellType = [wellType; true(np,1)];
            priIndex = [priIndex; (2+i)*ones(np,1)];
            gridSpacing = [gridSpacing; wellSpace];
            Pts= [Pts; wellPts];
            faultToCenter = [faultToCenter; zeros(size(wellPts,1),1)];
        else
            fracLine = linesToGrid{i};    
            [faultPts, fracSpace, CC, CR,CCid] = createFracGridPoints(fracLine,...
                                                                   faultGridSize,...
                                                                   opt.circleFactor);
            nl = size(faultPts,1)/2;
            if nl==0
                continue
            end

            newFaultType = lastFaultType+1:lastFaultType+nl;
            lastFaultType = newFaultType(end);
            faultType = [faultType; newFaultType';newFaultType'];        % 2*i-1rldecode([2*i-1; 2*i], [nl; nl])];
            wellType  = [wellType; false(2*nl,1)];
            priIndex = [priIndex; (2+i)*ones(2*nl,1)];
            gridSpacing = [gridSpacing;fracSpace];
            if isempty(faultPos)
                faultPos = size(Pts,1) + 1;
            end
            faultPos = [faultPos; size(Pts,1)+1+size(faultPts,1)];
            Pts = [Pts;faultPts];
            faultCenter = [faultCenter; CC]; % I think there is a bug here
            faultRadius = [faultRadius; CR]; % because wells and faults may
            faultToCenter = [faultToCenter;lastCCid+CCid]; % be interchanged
            lastCCid = faultToCenter(end)+1;

        end
    end
    

    %%
    %% Create reservoir grid
    dx = resGridSize;
    dy = dx;
    vx = 0:dx:pdims(1);
    vy = 0:dy:pdims(2);
    [X, Y] = meshgrid(vx, vy);

%    nx = numel(vx);
%    ny = numel(vy);
%    [ii, jj] = meshgrid(1:nx, 1:ny);
%     nedge = 1;
%     exterior = (ii <= nedge | ii > nx - nedge) | ...
%                (jj <= nedge | jj > ny - nedge);
% 
%     interior = ~exterior;
% 
%     X(interior) = X(interior);
%     Y(interior) = Y(interior);
     resPtsInit = [X(:), Y(:)];

    if any(wellType)
        res = {};
        varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', mlqtLevelSteps};
        for i = 1:size(resPtsInit,1)
            res = [res; mlqt(resPtsInit(i,:), Pts(wellType,:), dx, varArg{:})];
        end
        resPts = vec2mat([res{:,1}],2);
        resGridSize = 0.5*[res{:,2}]';
    else
        resPts = resPtsInit;
        resGridSize = repmat(max(dx,dy),size(resPts,1),1);
    end
    
    faultType = [faultType; zeros(size(resPts,1),1)];
    wellType = [wellType; false(size(resPts,1),1)];
    priIndex = [priIndex; max(priIndex) + ones(size(resPts, 1), 1)];
    gridSpacing = [gridSpacing; resGridSize]; %(1-10^-6)*min(dx,dy)*ones(size(resPts,1),1)];
    
    Pts = [Pts;resPts];
    
    % Remove Conflic Points
    [Pts, removed, wellType] = removeConflictPoints(Pts, gridSpacing, priIndex, wellType);
        
    faultType = faultType(~removed);
    wellType = wellType(~removed);
    
    if opt.fullFaultEdge    
        % Create a index map of the points
        faultIndex = (1:size(Pts,1))';
        remove = zeros(size(faultCenter,1),1);
        for i = 1:size(faultPos,1)-1 %iterate over all faults
            %Find removed fault points from current fault
            faultPts = faultIndex(faultPos(i):faultPos(i+1)-1);
            isRemovedFaults = false(size(Pts,1),1);
            isRemovedFaults(faultPts) = removed(faultPts);
            % Find the first and last circle in curren fault
            Cfrom  = faultToCenter(faultPos(i));
            Cto = faultToCenter(faultPos(i+1)-1)+1;
            % Find which circles has a removed fault point
            toRemove1 = false(size(faultCenter,1),1);
            toRemove2 = toRemove1;
            circleToRemove1 = faultToCenter(faultIndex(isRemovedFaults));
            circleToRemove1 = unique(circleToRemove1);
            % Remove circle if a fault ptn on both sides is removed
            circleToRemove2 = circleToRemove1 + 1;
            circleToRemove1 = [circleToRemove1; Cto];
            circleToRemove2 = [Cfrom; circleToRemove2];
            toRemove1(circleToRemove1) = true;
            toRemove2(circleToRemove2) = true;
            toRemove = and(toRemove1, toRemove2);
            remove = remove + toRemove;
        end
        remove = logical(remove);
        faultCenter = faultCenter(~remove,:);
        faultRadius = faultRadius(~remove);
        [Pts, removed] = enforceSufficientFaultCondition(Pts, faultType, faultCenter, faultRadius);
        faultType = faultType(~removed);
        wellType = wellType(~removed);
    end
    %% Create Grid
    Tri = delaunayTriangulation(Pts);
    G = triangleGrid(Pts, Tri.ConnectivityList);
    G = pebi(G);

    %wellType = wellType + cellsContPts(G, wellPts(removed(1:size(wellPts,1)),:));
    % label fault faces.
    N = G.faces.neighbors + 1;
    faultType = [0; faultType];
    ft1 = faultType(N(:,1));
    ft2 = faultType(N(:,2));
    G.faces.isFault = logical(ft1 == ft2 & ft1 > 0 & ft2 >0);
    
    %Label well cells
    G.cells.isWell = logical(wellType);

    varargout{1} = G;
    if nargout > 1
        varargout{2} = indicator;
    end
 
%% Plotting for debugging.
%         figure()
%     hold on
%     plot(Pts(logical(wellType),1),Pts(logical(wellType),2),'.')
%        %       hold on
%        plot(Pts(:,1),Pts(:,2),'r.')
% %       plot(fracPts(:,1), fracPts(:,2),'r.')
%       plotGrid(G,'facecolor','none')
%     
%     n = 50;
%     theta = (linspace(0,2*pi,n))';
%     for i = 1:size(Pts,1)
%         if wellType(i)
%         x = Pts(i,1) + gridSpacing(i)*cos(theta);
%         y = Pts(i,2) + gridSpacing(i)*sin(theta);
%         plot(x,y);
%         end
%     end
end


function [Pts, removed, wellType] = removeConflictPoints(Pts, gridSpacing, ...
                                                         priIndex, wellType)
    gridSpacing = gridSpacing*(1-1e-4); % To avoid floating point errors
    
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
        [~, Ii ] = sort(priIndex(Ic), 'descend');

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
