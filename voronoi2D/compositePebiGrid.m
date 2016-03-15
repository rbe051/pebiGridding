function G = compositePebiGrid(resGridSize, pdims, varargin)
    % Creates a PEBI grid adapting to faults and well traces.
    %
    % Argumets:
    %   resGridSize         Size of the reservoir grid cells
    %   pdims               [xmax, ymax], array with the size of the square
    %                       to be gridded.
    %
    % Varargin:
    %   wellLines           A struct of arrays. Each array is the 
    %                       coordinates of a well trace. If the
    %                       array only contains one coordinate, the well is
    %                       treated as a point well.
    %   wellGridFactor      The relative grid size of the well grid cells
    %                       compared to reservoir grid cells
    %   mlqtMaxLevel        Number of refinement steps to be used towards
    %                       wells
    %   mlqtLevelSteps      Array of size mlqtMaxLevel which specify the
    %                       radius of each refinement level
    %   faultLines          A struct of arrays. Each array is the
    %                       coordinates of a fault trace. 
    %   faultGridFactor     The relative grid size of the fault grid cells
    %                       compared to the reservoir grid cells
    %   circleFactor        Set the relative radius of the circles used to
    %                       create the fault cells.
    %   priOrder            Array of length number of wells + number of 
    %                       faults. Sets the priority of well and fault 
    %                       traces. First element set the priority of the
    %                       first well, last element set the priority of
    %                       the last fault.
    %   fullFaultEdge       Set to true if you wish to guarantee the faults
    %                       to be traced by edges in the PEBI grid
    %
    % Returns:
    %   G                   A mrst grid structure. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %% Set options
    opt = struct('wellLines',       {{}}, ...
                 'wellGridFactor',  0.5,  ...
                 'mlqtMaxLevel',    0,    ...
                 'mlqtLevelSteps',  -1,   ...
                 'faultLines',      {{}}, ...
                 'faultGridFactor', 0.5,  ...
                 'circleFactor',    0.6,  ...
                 'fullFaultEdge',   0,    ...
                 'priOrder',        []);         
    
    opt = merge_options(opt, varargin{:});
    circleFactor = opt.circleFactor;
    
    % Set grid sizes
    wellGridSize   = resGridSize*opt.wellGridFactor;
    faultGridSize  = resGridSize*opt.faultGridFactor;
    mlqtMaxLevel   = opt.mlqtMaxLevel;
    mlqtLevelSteps = opt.mlqtLevelSteps;
    
    % Load faults and Wells
    faultLines  = opt.faultLines;
    wellLines   = opt.wellLines;
    nFault      = numel(faultLines);
    nWell       = numel(wellLines);
    linesToGrid = [wellLines, faultLines{:}];
    isWell      = [true(nWell,1); false(nFault,1)];
    
    % Set priority index
    priOrder = opt.priOrder;
    if isempty(priOrder)
        priOrder = 1:nFault+nWell;
    else
        assert(numel(unique(priOrder)) == nFault + nWell);
    end
    % Sort faults and wells in priority
    linesToGrid = linesToGrid(priOrder);
    isWell = isWell(priOrder);
    
    
    %% Test input
    assert(resGridSize>0);
    assert(numel(pdims)==2);
    assert(all(pdims>0 ));
    assert(wellGridSize>0);
    assert(mlqtMaxLevel>=0);
    assert(faultGridSize>0);
    assert(0.5<circleFactor && circleFactor<1);
    

    %% Initialize variables.
    faultType   = [];           % To keep track of fault points
    wellType    = logical([]);  % To keep track of well points 
    priIndex    = [];           % Priority of each point
    gridSpacing = [];           % Allowed grid spacing for each point

    Pts           = [];         % Points
    faultCenter   = [];         % Center of circle used to create fault pts
    faultRadius   = [];         % Radius of the circle
    faultToCenter = [];         % Map from a fault to the circle center
    faultPos      = [];         % Start and end index of a fault in Pts


    %% Create fault and well points
    lastFaultType = 0;
    lastCCid = 0;
    for i = 1:nFault + nWell  % create well points
        if isWell(i)
            wellLine      = linesToGrid{i};
            [wellPts,wellSpace] = createWellGridPoints(wellLine, wellGridSize);
            
            np            = size(wellPts,1);
            faultType     = [faultType; zeros(np,1)];
            wellType      = [wellType; true(np,1)];
            priIndex      = [priIndex; (2+i)*ones(np,1)];
            gridSpacing   = [gridSpacing; wellSpace];
            Pts           = [Pts; wellPts];
            faultToCenter = [faultToCenter; zeros(size(wellPts,1),1)];
        end
    end
    for i = 1:nFault + nWell % create fault points
        if ~isWell(i)
            fracLine      = linesToGrid{i};    
            [faultPts, fracSpace, CC, CR, CCid] =                       ...
                                    createFaultGridPoints(fracLine,     ... 
                                                          faultGridSize,...
                                                          circleFactor);
            nl = size(faultPts,1)/2;
            if nl==0 % No fault points created
                continue
            end

            newFaultType  = lastFaultType+1:lastFaultType+nl;
            lastFaultType = newFaultType(end);
            faultType     = [faultType; newFaultType';newFaultType'];
            wellType      = [wellType; false(2*nl,1)];
            priIndex      = [priIndex; (2+i)*ones(2*nl,1)];
            gridSpacing   = [gridSpacing;fracSpace];
            if isempty(faultPos)
                faultPos = size(Pts,1) + 1;
            end
            faultPos      = [faultPos; size(Pts,1)+1+size(faultPts,1)];
            Pts           = [Pts;faultPts];
            faultCenter   = [faultCenter; CC];
            faultRadius   = [faultRadius; CR]; 
            faultToCenter = [faultToCenter;lastCCid+CCid];
            lastCCid      = faultToCenter(end)+1;

        end
    end
    

    %% Create reservoir grid
    dx = pdims(1)/ceil(pdims(1)/resGridSize);
    dy = pdims(2)/ceil(pdims(2)/resGridSize);
    vx = 0:dx:pdims(1);
    vy = 0:dy:pdims(2);

    [X, Y] = meshgrid(vx, vy);

    resPtsInit = [X(:), Y(:)]; % Reservoir grid points before refinement

    
    % Refine reservoir grid
    if any(wellType)
        res = {};
        varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', mlqtLevelSteps};
        for i = 1:size(resPtsInit,1)
            res = [res; mlqt(resPtsInit(i,:), Pts(wellType,:), resGridSize, varArg{:})];
        end
        resPts      = vec2mat([res{:,1}],2);
        resGridSize = 0.5*[res{:,2}]';
    else
        resPts      = resPtsInit;
        resGridSize = repmat(0.5*min(dx,dy),size(resPts,1),1);
    end
    
    faultType   = [faultType; zeros(size(resPts,1),1)];
    wellType    = [wellType; false(size(resPts,1),1)];
    priIndex    = [priIndex; max(priIndex) + ones(size(resPts, 1), 1)];
    gridSpacing = [gridSpacing; resGridSize];
    
    Pts = [Pts;resPts];
    
    %% Remove Conflic Points
    [Pts, wellType, removed] = removeConflictPoints(Pts, gridSpacing, priIndex, wellType);
        
    faultType = faultType(~removed);
    
    if opt.fullFaultEdge
        [faultCenter, faultRadius] = removeFaultCircles(faultCenter, ...
                                                        faultRadius,...
                                                        faultPos,...
                                                        faultToCenter,...
                                                        removed);
        [Pts, removed, ~, ~] = enforceSufficientFaultCondition(Pts, ...
                                                               faultType,...
                                                               faultCenter,...
                                                               faultRadius);
        % If you whish to use the updated faultCenter and faultRadius later 
        % you can retrive them by replacing the ~ outputs.
        
        faultType = faultType(~removed);
        wellType = wellType(~removed);
    end
    
    
    %% Create Grid
    Tri = delaunayTriangulation(Pts);
    G = triangleGrid(Pts, Tri.ConnectivityList);
    G = pebi(G);

    % label fault faces.
    N = G.faces.neighbors + 1;
    faultType = [0; faultType];
    ft1 = faultType(N(:,1));
    ft2 = faultType(N(:,2));
    G.faces.tag = logical(ft1 == ft2 & ft1 > 0 & ft2 >0);
    
    %Label well cells
    G.cells.tag = logical(wellType);

end


function [CC, CR] = removeFaultCircles(CC, CR,faultPos,faultToCenter, remFaultPts)
    assert(size(CC,1)==size(CR,1));    
    % Create a index map of the points
    faultIndex = (1:size(faultToCenter,1))';
    remCirc = zeros(size(CC,1),1);
    for i = 1:size(faultPos,1)-1 %iterate over all faults
        % Find removed fault points from current fault
        faultPts = faultIndex(faultPos(i):faultPos(i+1)-1);
        isRemovedFaults = false(size(faultToCenter,1),1);
        isRemovedFaults(faultPts) = remFaultPts(faultPts);
        % Find the first and last circle in curren fault
        Cfrom  = faultToCenter(faultPos(i));
        Cto = faultToCenter(faultPos(i+1)-1)+1;
        % Find which circles has a removed fault point
        toRemove1 = false(size(CC,1),1);
        toRemove2 = toRemove1;
        circleToRemove1 = faultToCenter(faultIndex(isRemovedFaults));
        circleToRemove1 = unique(circleToRemove1);
        % Remove circle if a fault pnt on both sides is removed
        circleToRemove2 = circleToRemove1 + 1;
        circleToRemove1 = [circleToRemove1; Cto];
        circleToRemove2 = [Cfrom; circleToRemove2];
        toRemove1(circleToRemove1) = true;
        toRemove2(circleToRemove2) = true;
        toRemove = and(toRemove1, toRemove2);
        remCirc = remCirc + toRemove;
    end
    remCirc = logical(remCirc);
    CC = CC(~remCirc,:);
    CR = CR(~remCirc);    

end

function [Pts, removed, CC, CR] = ...
            enforceSufficientFaultCondition(Pts, isFault, CC, CR)

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
