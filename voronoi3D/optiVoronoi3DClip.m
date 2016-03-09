function [G,optPts,f,g] = optiVoronoi3DClip(pts, dt, varargin)
    % Creates an optimal voronoi voronoi diagram. By optimal means that all
    % voronoi seeds coincide with the cell mass centers.
    %
    % Arguments:
    %   pts         n x 3 array of initial of seeds points
    %   bndr        nb x 3 array of points defining the convext boundary
    % Varagin:
    %   density     density function in the domain
    %   storedVec   number of stored vectors used to approximate the
    %               hessian in the lbfgs algorithm
    %   tol         convergence tolerance
    %   maxIt       maximum number of iteration in the lbfgs algorithm
    %   minStep     if step length is smaler than minStep, the lbfgs
    %               algorithm stops
    % Returns:
    %   G           mrst grid structure
    %   optPts      optimal grid seeds
    %   f           function value at each step
    %   g           gradient value at each step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Runar Lie Berge                                                  2016
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    assert (size(pts, 2) == 3, ...
          ['Function ''%s'' is only supported in three ', ...
           'space dimensions.'], mfilename);
 

    opt = struct('density',   @(x) ones(size(x,1),1),...
                 'fault', struct,                    ...
                 'storedVec', 5,                     ...
                 'tol',       1e-10,                 ...         
                 'maxIt',     1000,                  ...
                 'minStep',   10*eps);
    opt = merge_options(opt, varargin{:});
    
    F = @(pts) objectiveFunc(pts, dt, opt.fault, opt.density);

    pts = reshape(pts',[],1);

    [optPts, f, g] = lbfgs(pts, F, dt, 'storedVec',opt.storedVec, ...
                                       'maxIt',    opt.maxIt,...
                                       'tol',      opt.tol);

    optPts = reshape(optPts,3,[])';
    G = restrictedVoronoiDiagram(optPts, dt);
end


function [f, g] = objectiveFunc(pts, bndr, fault, rho)
    pts = reshape(pts,3,[])';

    G = restrictedVoronoiDiagram(pts, bndr);
    G = computeGeometry(G);
    dt = delaunayTriangulation(pts);
    E = dt.edges();
    [V,C,symV] = clipGrid(dt,fault);
    [n,x0] = normalPlanes(V,C); %OBS! This is simplified and subject to change.
    [V,C,cells, symV] = VOuter(V,C,symV,pts,G,dt,n,x0);
    [T, triPos, vol] = triang(V,C,cells);

    gv = volumeGrad(V,T,triPos,cells,symV,pts,E,n, fault);
    
    intFun = @(x,i) sum(repmat(rho(x),1,3).*(x-repmat(pts(i,:),size(x,1),1)).^2,2);
    
    
    f = sum(polyhedronInt(G,1:G.cells.num, intFun));

    massFun = @(x,i) rho(x);
    masses = polyhedronInt(G,1:G.cells.num,massFun);
    g = reshape((2*repmat(masses,1,3).*(pts - G.cells.centroids))',[],1);

end


function [g] = volumeGrad(V,T,triPos,cells,symV,pts,E,n,fault)
  g = zeros(size(pts));
  for i = 1:numel(triPos)-1
    tri = T(triPos(i):triPos(i+1)-1,:);
    U = V(tri,:)-repmat(pts(i,:),size(V,1),1);
    for j = 1:size(tri,1)
      dvds = zeros(1,3);
      for k = 1:size(tri,2)
        Uact = U(tri(j,1:3~=k),:);
        dvdc = 1/6*cross(Uact(1,:),Uact(2,:),1);
        dvds = dvds - dvdc;
        switch sum(symV{tri(j,k)}>0)
          case 0
            continue
          case 1
            bis = symV{tri(j,k)}(symV{tri(j,k)}>0);
            %f =  symV{tir(i)}(symV{tri(i)}<0);
            n2 = [1,1,0]; %OBS This MUST be made general
            s0 = cells(i);          
            s2 = E(bis,E(bis,:)~=s0);

            A = [(pts(s2,:)-pts(s0));n;n2];
            B = [(V(tri(j,k),:)-pts(s0,:)),(pts(s2,:)-V(tri(j,k),:)) ;...
                 zeros(2,6)];
          case 2
            bis = symV{tri(j,k)}(symV{tri(j,k)}>0);
            s0= cells(i);
            e = unique(E(bis,:));
            s = e(e~=s0);
            s2 = s(1); s3 = s(2);
            A = [(pts(s2,:)-pts(s0));(pts(s3,:)-pts(s0));n];
            B = [V(tri(j,k),:)-pts(s0,:),(pts(s2,:)-V(tri(j,k),:)),zeros(1,3);...
                 V(tri(j,k),:)-pts(s0,:), zeros(1,3), pts(s3,:)-V(tri(j,k),:);...
                 zeros(1,9)];
            [e,I] = sort(e);
            gradC = A\B;
            gradC = [zeros(3,e(1)-1)     ,gradC(:,3*I(1)-2:3*I(1)),...
                     zeros(3,e(2)-e(1)-1),gradC(:,3*I(2)-2:3*I(2)),...
                     zeros(3,size(pts,1)-e(2))];
          case 3
            bis = symV{tri(j,k)};
            s0 = cells(i);

            e = unique(E(bis,:));
            s = e(e~=s0);
            s2 = s(1); s3 = s(2); s4 = s(3);
            A = [(pts(s2,:)-pts(s0));pts(s3,:)-pts(s0);pts(s4,:)-pts(s0)];
            B = [V(tri(j,k),:)-pts(s0,:),(pts(s2,:)-V(tri(j,k),:)),zeros(1,6);...
                 V(tri(j,k),:)-pts(s0,:), zeros(1,3), pts(s3,:)-V(tri(j,k),:),zeros(1,3);...
                 V(tri(j,k),:)-pts(s0,:), zeros(1,6), pts(s4,:)-V(tri(j,k),:)];
            [e,I] = sort(e);
            gradC = A\B;
            gradC = [zeros(3,e(1)-1)     ,gradC(:,3*I(1)-2:3*I(1)),...
                     zeros(3,e(2)-e(1)-1),gradC(:,3*I(2)-2:3*I(2)),...
                     zeros(3,e(3)-e(2)-1),gradC(:,3*I(3)-2:3*I(3)),...
                     zeros(3,size(pts,1)-e(3))];
          otherwise
            warning('this should not happen!')
        end
        
        g = g + dvdc*gradC;
      end
      g = g + [zeros(1,cells(i)),dvds,zeros(1,size(pts,3)-ce)
    end
  end

end


function [T,triPos, vol] = triang(V,C,cells)
T = [];
triPos = zeros(numel(cells)+1,1);
triPos(1)=1;
vol = 0;
for i = 1:numel(cells)
    c = C{cells(i)};
    [hull,volAdd] = convhull(V(c,:));
    T = [T;c(hull)]; 
    triPos(i+1) = triPos(i) + size(hull,1);  
    vol = vol+volAdd;
end
end


function [n,x0] = normalPlanes(V,C)

x0 = cellfun(@(c) mean(V(c,:)),C,'uniformOutput',false);
x0 = cell2mat(x0);

vert = cell2mat(C');
n = cross(V(vert(2),:)-V(vert(1),:),V(vert(3),:)-V(vert(1),:));
n = n/norm(n,2);
% for i =1:numel(C)  % This is to approximate a mean plane. For now I
% assume that the fault is a plane
%     p = V(C{i},:)- x0(i,:);
%     [~,~,N] = svd(p);
%     n(i) = N(:,end)';
% end

end


function [V,C,cells, symV] = VOuter(V,C,symV,s,G,dt,n,x0)
cells = [];
for i=1:numel(C)
    if ~isempty(C{i})
        faces = G.cells.faces(G.cells.facePos(i):G.cells.facePos(i+1)-1);
        nodes = arrayfun(@(l,r) G.faces.nodes(l:r), ...
                                G.faces.nodePos(faces),...
                                G.faces.nodePos(faces+1)-1, 'un', 0)';
        nodes = unique(cell2mat(nodes'));
        outer = sign((bsxfun(@minus,G.nodes.coords(nodes,:),x0(i,:)))*n') ...
              ~=sign((s(i,:)-x0(i,:))*n');
        C{i} = [C{i},size(V,1)+1:size(V,1)+sum(outer)];
        V = [V; G.nodes.coords(nodes(outer),:)];        
        cells = [cells;i];
        symV = [symV;findBisector(G,dt.edges,nodes(outer),i)];
    end
end
end



function [b] = findBisector(G,dtEdges,nodes,c)

b = cell(numel(nodes),1);
for i = 1:numel(nodes)
   [~, face] = max(bsxfun(@gt,G.faces.nodePos,find(G.faces.nodes==nodes(i))'),[],1);
   face = face-1;
   cells = sort(G.faces.neighbors(face,:),2);
   [~,b{i}] = ismember(cells(any(cells==c,2),:), dtEdges,'rows');
end
    
end





































