function varargout = compositeGridPEBIdistmesh(resGridSize, pdims, varargin)
        
    opt = struct('wellLines', {{}}, ...
                 'wellGridFactor', -1, ...
                 'wellRefDist', -1,...
                 'faultGridFactor', -1, ...
                 'faultLines', {{}}, ...
                 'circleFactor', 0.6, ...
                 'priOrder', []);
             
    opt = merge_options(opt, varargin{:});

    %% Set options
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
    wellEps = opt.wellRefDist; % For grid refinement.
    if wellEps<0
        wellEps = 0.25/max(pdims);
    end
    
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

    fixedPts = [];
    %faultCenter = [];
    %faultRadius = [];
    %faultToCenter = zeros(size(wellPts,1),1);
    %faultPos = size(wellPts,1)+1;


    %% Place fault and well points
    lastFaultType = 0;
    %lastCCid = 0;
    for i = 1:nFault+nWell % From low priority to high
        if isWell(i)
            wellLine = linesToGrid{i};
            [wellPts,wellSpace] = createWellGridPoints(wellLine, wellGridSize);
            np = size(wellPts,1);
            faultType = [faultType; zeros(np,1)];
            wellType = [wellType; true(np,1)];
            priIndex =[priIndex; i*ones(2*np,1)];
            gridSpacing = [gridSpacing; wellSpace];
            fixedPts= [fixedPts; wellPts];
        end
    end
    %% create distance function
    if nWell>0
        h = @(x) min((ones(size(x,1),1)/wellGridFactor), ...
                     min(1.2*exp(pdist2(x,fixedPts(wellType,:))/wellEps),[],2));
    else
        h = @(p) huniform(p)/wellGridFactor;
    end
    hfault = @(x) min((ones(size(x,1),1)*faultGridFactor/wellGridFactor), ...
                     min(1.2*exp(pdist2(x,fixedPts(wellType,:))/wellEps),[],2));
    for i = 1:nFault + nWell

       if ~isWell(i)
           fracLine = linesToGrid{i};    
           [faultPts, fracSpace,~,~,~] = createFracGridPoints(fracLine,...
                                                              wellGridSize,...
                                                              opt.circleFactor,...
                                                              'distFunc', hfault);
           nl = size(faultPts,1)/2;
           if nl==0
               continue
           end

           newFaultType = lastFaultType+1:lastFaultType+nl;
           lastFaultType = newFaultType(end);
           faultType = [faultType; newFaultType';newFaultType'];        % 2*i-1rldecode([2*i-1; 2*i], [nl; nl])];
           wellType  = [wellType; false(2*nl,1)];
           priIndex = [priIndex; i*ones(2*nl,1)]; 
           gridSpacing = [gridSpacing;fracSpace];
           fixedPts = [fixedPts;faultPts];
           %faultCenter = [faultCenter; CC];
           %faultRadius = [faultRadius; CR];
           %faultToCenter = [faultToCenter;lastCCid+CCid]; 
           %lastCCid = faultToCenter(end)+1;
           %faultPos = [faultPos; faultPos(end)+size(newFracPts,1)]; 
       end
    end
    

    %% Remove fault and well conflic points
    if size(fixedPts,1)>1
        [fixedPts,removed,wellType]=removeConflictPoints(fixedPts, ...
                                                         gridSpacing,...
                                                         priIndex, ...
                                                         wellType);        
        wellType = wellType(~removed);
        faultType = faultType(~removed);
        gridSpacing = gridSpacing(~removed);
        priIndex = priIndex(~removed);
    end
    

    %%
    %% Create grid
    % set dist function
    x = pdims;
    fd = @(p) drectangle(p,0,x(1),0,x(2));
    % Set fixed points
    corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];
    fixedPts = [fixedPts; corners]; %Add these to fault and well type later
    
    [Pts,t,sort] = distmesh2d(fd, h, wellGridSize,[0,0;x(1),x(2)], fixedPts);
    nNewPts = size(Pts,1) - size(faultType,1);
    
    faultType = [faultType;zeros(nNewPts,1)];
    wellType = [wellType;zeros(nNewPts,1)];
    gridSpacing = [gridSpacing;zeros(nNewPts,1)];
    priIndex = [priIndex; max(priIndex) + ones(nNewPts,1)];
    
    faultType = faultType(sort);
    wellType = wellType(sort);
    gridSpacing = gridSpacing(sort);
    priIndex = priIndex(sort);
    %gridSpacing(logical(wellType)) = gridSpacing(logical(wellType));

    %% Remove new conflict  points
    [Pts,removed,wellType]=removeConflictPoints(Pts, ...
                                                gridSpacing,...
                                                priIndex, ...
                                                wellType);        
    wellType = wellType(~removed);
    faultType = faultType(~removed);
    %gridSpacing = gridSpacing(~removed); % Should be added of used later
    %priIndex = priIndex(~removed);
    t = delaunay(Pts);
    pmid=(Pts(t(:,1),:)+Pts(t(:,2),:)+Pts(t(:,3),:))/3;    % Compute centroids
    t=t(feval(fd,pmid)<-0.001*wellGridFactor,:);  % Keep interior triangles

    [Pts,t, ~, sort]=fixmesh(Pts,t);
    faultType = faultType(sort);
    wellType = wellType(sort);
    %gridSpacing = gridSpacing(sort);
    %priIndex = priIndex(sort);
    %%
    G = triangleGrid(Pts, t);
    G = pebi2(G);
    

    %wellType = wellType + cellsContPts(G, wellPts(removed(1:size(wellPts,1)),:));
    %label fault faces.
    N = G.faces.neighbors + 1;
    faultType = [0; faultType];
    ft1 = faultType(N(:,1));
    ft2 = faultType(N(:,2));
    G.faces.isFault = logical(ft1==ft2 & ft1 > 0 & ft2 > 0);
    
    %Label well cells
    G.cells.isWell= logical(wellType);

    varargout{1} = G;
    if nargout > 1
        varargout{2} = indicator;
    end
    
    
    %% Plotting for debugging.
%     figure()
%      hold on
% %     plot(Pts(logical(wellType),1),Pts(logical(wellType),2),'.')
% %        figure()
% %     hold on
% %        %       hold on
% %       plot(Pts(:,1),Pts(:,2),'r.')
%     plot(fixedPts(:,1), fixedPts(:,2),'r.')
%     plotGrid(G,'facecolor','none')
%     
%     n = 50;
%     theta = (linspace(0,2*pi,n))';
%     for i = 1:size(fixedPts,1)-size(corners,1)
%         x =  fixedPts(i,1) + gridSpacing(i)*cos(theta);
%         y = fixedPts(i,2) + gridSpacing(i)*sin(theta);
%         plot(x,y);     
%     end
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