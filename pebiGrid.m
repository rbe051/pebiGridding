function G = pebiGrid(resGridSize, pdims, varargin)
    % Creates a PEBI grid adapting to faults and well traces.
    %
    % Argumets:
    %   resGridSize         Size of the reservoir grid cells
    %   pdims               [xmax, ymax], array with the size of the square
    %                       to be gridded.
    %
    % Varargin:
    %   wellLines           A struct of arrays. Each array is the 
    %                       coordinates of a well trace. If an
    %                       array only contains one coordinate, the well is
    %                       treated as a point well.
    %   wellGridFactor      The relative grid size of the well grid cells
    %                       compared to reservoir grid cells
    %   wellRefinement      Logical scalar set to true if refinement
    %                       towards wells is wanted
    %   epsilon             Scale for grid refinement.
    %   faultLines          A struct of arrays. Each array is the
    %                       coordinates of a fault trace. 
    %   faultGridFactor     The relative grid size of the fault grid cells
    %                       compared to the reservoir grid cells
    %   circleFactor        The relative radius of the circles used to
    %                       create the fault cells.
    %   priOrder            Array of length = number of wells + number of 
    %                       faults. Sets the priority of well and fault 
    %                       traces. First element set the priority of the
    %                       first well, last element set the priority of
    %                       the last fault.
    %   fullFaultEdge       Set to true if you wish to guarantee the faults
    %                       to be traced by edges in the PEBI grid
    %
    % Returns:
    %   G                   A mrst grid structure. 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Runar Lie Berge (runarlb@stud.ntnu.no)
%% January 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

    %% Set options
    opt = struct('wellLines',       {{}}, ...
                 'wellGridFactor',  0.5, ...
                 'wellRefinement',  false, ...
                 'epsilon',         -1,...
                 'faultLines',      {{}}, ...
                 'faultGridFactor', 0.5, ...
                 'circleFactor',    0.6, ...
                 'priOrder',        []);  
             
    opt          = merge_options(opt, varargin{:});
    circleFactor = opt.circleFactor;
	wellRef      =  opt.wellRefinement;
    wellEps      = opt.epsilon; % For grid refinement.
    % Set grid sizes
    wellGridFactor  = opt.wellGridFactor;
    faultGridFactor = opt.faultGridFactor;
    wellGridSize    = resGridSize*wellGridFactor;
    faultGridSize   = resGridSize*faultGridFactor;

    if wellEps<0
        wellEps = 0.25/max(pdims);
    end
    
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
        assert(numel(priOrder) == nFault + nWell);
    end
    % Sort faults and wells in priority
    linesToGrid = linesToGrid(priOrder);
    isWell      = isWell(priOrder);
    
    %% Test input
    assert(resGridSize>0);
    assert(numel(pdims)==2);
    assert(all(pdims>0 ));
    assert(wellGridSize>0);
    assert(faultGridSize>0);
    assert(0.5<circleFactor && circleFactor<1);
    
    %% Initialize variables.
    faultType   = [];
    wellType    = logical([]);
    priIndex    = [];
    gridSpacing = [];
    fixedPts = [];

    %% Place well points
    lastFaultType = 0;
    %lastCCid = 0;
    for i = 1:nFault+nWell % From low priority to high
        if isWell(i)
            wellLine = linesToGrid{i};
            [wellPts,wellSpace] = createWellGridPoints(wellLine, wellGridSize);
            
            np          = size(wellPts,1);
            faultType   = [faultType; zeros(np,1)];
            wellType    = [wellType; true(np,1)];
            priIndex    = [priIndex; i*ones(np,1)];
            gridSpacing = [gridSpacing; wellSpace];
            fixedPts    = [fixedPts; wellPts];
        end
    end
    %% create distance function
    if wellRef && nWell>0
        hres   = @(x) min((ones(size(x,1),1)/wellGridFactor), ...
                      min(1.2*exp(pdist2(x,fixedPts(wellType,:))/wellEps),[],2));
        hfault = @(x) wellGridSize*faultGridFactor*hres(x);
    else
    hres   = @(p) constFunc(p)/wellGridFactor;
	hfault = @(p) constFunc(p)*faultGridSize;
    end
    
    %% Place fault points
    for i = 1:nFault + nWell
       if ~isWell(i)
           fracLine = linesToGrid{i};    
           [faultPts, fracSpace,~,~,~] = createFaultGridPoints(fracLine,...
                                                               faultGridSize,...
                                                               circleFactor,...
                                                              'distFunc', hfault);
           nl = size(faultPts,1)/2;
           if nl==0
               continue
           end

           newFaultType  = lastFaultType+1:lastFaultType+nl;
           lastFaultType = newFaultType(end);
           faultType     = [faultType; newFaultType';newFaultType'];
           wellType      = [wellType; false(2*nl,1)];
           priIndex      = [priIndex; i*ones(2*nl,1)]; 
           gridSpacing   = [gridSpacing;fracSpace];
           fixedPts      = [fixedPts;faultPts];
       end
    end
    

    %% Remove fault and well conflic points
    if size(fixedPts,1)>1
        [fixedPts, wellType, removed]=removeConflictPoints(fixedPts, ...
                                                           gridSpacing,...
                                                           priIndex, ...
                                                           wellType);        
        faultType   = faultType(~removed);
        gridSpacing = gridSpacing(~removed);
        priIndex    = priIndex(~removed);
    end
    


    %% Create Reservoir grid points
    % set dist function
    x = pdims;
    rectangle = [0,0; x(1),x(1)];   
    fd = @(p) drectangle(p, 0, x(1), 0, x(2));
    % Set fixed points
    corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];


    fixedPts = [fixedPts; corners];
    [Pts,~,sorting] = distmesh2d(fd, hres, wellGridSize, rectangle, fixedPts);

    nNewPts     = size(Pts,1) - size(faultType,1);
    faultType   = [faultType;zeros(nNewPts,1)];
    wellType    = [wellType;zeros(nNewPts,1)];
    gridSpacing = [gridSpacing;zeros(nNewPts,1)];
    priIndex    = [priIndex; max(priIndex) + ones(nNewPts,1)];
    
    Pts         = Pts(sorting,:);
    faultType   = faultType(sorting);
    wellType    = wellType(sorting);
    gridSpacing = gridSpacing(sorting);
    priIndex    = priIndex(sorting);
    
    %% Remove new conflict  points
    [Pts, wellType, removed,]=removeConflictPoints(Pts, ...
                                                   gridSpacing,...
                                                   priIndex, ...
                                                   wellType);        
    faultType = faultType(~removed);
    %gridSpacing = gridSpacing(~removed); % Should be added if used later
    %priIndex = priIndex(~removed);
    
    %% Create grid
    t    = delaunay(Pts);
    % Fix boundary
    pmid = (Pts(t(:,1),:)+Pts(t(:,2),:)+Pts(t(:,3),:))/3;% Compute centroids
    t    = t(feval(fd,pmid)<-0.001*wellGridFactor,:);    % Keep interior triangles

    [Pts,t, ~, sorting] = fixmesh(Pts,t);
    faultType           = faultType(sorting);
    wellType            = wellType(sorting);
    %gridSpacing = gridSpacing(sort); % Should be added if used later
    %priIndex = priIndex(sort);
    
    G = triangleGrid(Pts, t);
    G = pebi(G);
    G = computeGeometry(G);
    
    %label fault faces.
    N               = G.faces.neighbors + 1;
    faultType       = [0; faultType];
    ft1             = faultType(N(:,1));
    ft2             = faultType(N(:,2));
    G.faces.tag = logical(ft1==ft2 & ft1 > 0 & ft2 > 0);
    
    %Label well cells
    G.cells.tag= logical(wellType);    
    
end

