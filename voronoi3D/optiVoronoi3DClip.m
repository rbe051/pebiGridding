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
                 'fault',     struct,                ...
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
    bndr2D.ConnectivityList = bndr.freeBoundary;
    bndr2D.Points = bndr.Points;
    [G,Gt,symG] = restrictedVoronoiDiagram(pts, bndr);
    G = computeGeometry(G);
    dt = delaunayTriangulation(pts);
    E = dt.edges();
    [nF,x0F]     = normalOfThreePts(reshape(fault.Points(fault.ConnectivityList',:)',9,[])');
    [nBdr,x0bdr] = normalOfThreePts(reshape(bndr.Points(bndr.freeBoundary',:)',9,[])');
    
    [Vf,Cf,symF] = clipGrid(dt,fault,'bisectN',nBdr,'bisectX0',x0bdr);

    %[n,x0] = normalPlanes(Vf,Cf); %OBS! This is simplified and subject to change.
    n = [nF;nBdr]; x0 = [x0F;x0bdr];
    
    cells = find(~cellfun(@isempty,Cf));
    [Vb,Cb, symB] = VOuter(cells,symG,pts,Gt,G,dt,n(1,:),x0(1,:));

    
    [V,C,symV,cutT] = mergeNodes(Vf,Cf,symF,fault,Vb,Cb,symB,bndr2D);

    [T, triPos, vol] = triang(V,C,cells);
    
    gv = volumeGrad(G,V,T,triPos,cells,symV,pts,E,n, cutT)';
    
    intFun = @(x,i) sum(repmat(rho(x),1,3).*(x-repmat(pts(i,:),size(x,1),1)).^2,2);
    
    
    f = sum(polyhedronInt(G,1:G.cells.num, intFun));
    f = 0*f + vol;

    massFun = @(x,i) rho(x);
    masses = polyhedronInt(G,1:G.cells.num,massFun);
    g = reshape((2*repmat(masses,1,3).*(pts - G.cells.centroids))',[],1);
    g = 0*g + 1*gv;
end




function [edg] = isOnEdge(dt,f, V)
  vert = dt.ConnectivityList(f,:);
  edg = false(size(vert,2),1);
  for i = 1:size(edg,1)
    edg(i) = isColinear([dt.Points(vert(i),:);V;dt.Points(vert(1+mod(i,size(edg,1))),:)]);
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


function [V,C,symV] = VOuter(c,symG,s,Gt,G,dt,n,x0)

symV = {};
C = cell(G.cells.num,1);
V = [];
for i=1:numel(c)
  faces = G.cells.faces(G.cells.facePos(c(i)):G.cells.facePos(c(i)+1)-1);
  nodes = arrayfun(@(l,r) G.faces.nodes(l:r), ...
                          G.faces.nodePos(faces),...
                          G.faces.nodePos(faces+1)-1, 'un', 0)';
  nodes = unique(cell2mat(nodes'));
  outer = sign((bsxfun(@minus,G.nodes.coords(nodes,:),x0))*n') ...
        ~=sign((s(c(i),:)-x0)*n');
  C{c(i)} = size(V,1)+1:size(V,1)+sum(outer);
  V = [V; G.nodes.coords(nodes(outer),:)];        
  symV = [symV;findBisector(G,Gt,symG,dt.edges,nodes(outer),c(i))];
end
    [V,IA,IC] = uniquetol(V,50*eps,'byRows',true);
    symV = symV(IA);
    C = cellfun(@(c) unique(IC(c))', C,'UniformOutput',false);
end


function [V,C,symV,dt] = mergeNodes(Vf,Cf,symF,F,Vb,Cb,symB,B)
nff = size(F.ConnectivityList,1);
ncf = numel(Cf);
V = [Vf;Vb];
C = cellfun(@(c1,c2) [c1,c2+size(Vf,1)], Cf, Cb,'un',false);
symV = [symF; cellfun(@(b) b - nff*(1-sign(b))/2, symB,'un',false)];
[V,IA,IC] = uniquetol(V,50*eps,'byRows',true);
symV = symV(IA);
C = cellfun(@(c) unique(IC(c))', C,'UniformOutput',false);

dt.Points = [F.Points;B.Points];
dt.ConnectivityList = [F.ConnectivityList; B.ConnectivityList + size(F.Points,1)];

end


function [b] = findBisector(G,Gt, symG,dtEdges,nodes,c)

b = cell(numel(nodes),1);
for i = 1:numel(nodes)
  if isempty(symG{nodes(i)})
   [~, face] = max(bsxfun(@gt,G.faces.nodePos,find(G.faces.nodes==nodes(i))'),[],1);
   face = face-1;
   cells = sort(G.faces.neighbors(face,:),2);
   [~,b{i}] = ismember(cells(any(cells==c,2),:), dtEdges,'rows');
  else
    b{i} = symG{nodes(i)};
  end
end
    
end




function [g] = volumeGrad(G,V,T,triPos,cells,symV,pts,E,nF,bndr)
  g = zeros(1,numel(pts));
  vol = 0;
  for i = 1:numel(triPos)-1
    tri = T(triPos(i):triPos(i+1)-1,:);
    U = V(tri',:)-repmat(pts(i,:), numel(tri),1);
    for j = 1:size(tri,1)
      dvds = zeros(1,3);
      vol =vol+ 1/6*dot(U(3*j-2,:), cross(U(3*j-1,:),U(3*j,:)));
      for k = 1:3
        switch k
          case 1
            Uact = U([3*j-1,3*j],:);
          case 2
            Uact = U([3*j,3*j-2],:);
          case 3 
            Uact = U([3*j-2,3*j-1],:); 
        end
        dvdc = 1/6*cross(Uact(1,:),Uact(2,:),2);
        dvds = dvds - dvdc;
        switch sum(symV{tri(j,k)}>0)
          case 0
            continue
          case 1
            bis = symV{tri(j,k)}(symV{tri(j,k)}>0);
            f =  symV{tri(j,k)}(symV{tri(j,k)}<0);
            assert(numel(f)==2);
            %n = matrixNormals(bndr,-f(1),V(tri(j,k),:));
            n  = nF(-f,:);
            
            s0 = cells(i);
            e = unique(E(bis,:));
            I = 1:numel(e);
            s0Id = e==s0;
            I = [I(s0Id), I(~s0Id)];
            s1 = e(~s0Id);

            A = [(pts(s1,:)-pts(s0,:));n];
            B = [(V(tri(j,k),:)-pts(s0,:)),(pts(s1,:)-V(tri(j,k),:)) ;...
                 zeros(2,6)];
            gradC = A\B;
            gradC = [zeros(3,3*(e(1)-1))     , gradC(:,3*I(1)-2:3*I(1)),...
                     zeros(3,3*(e(2)-e(1)-1)), gradC(:,3*I(2)-2:3*I(2)),...
                     zeros(3,numel(pts)-3*e(2))];
          case 2
            bis = symV{tri(j,k)}(symV{tri(j,k)}>0);
            f =  symV{tri(j,k)}(symV{tri(j,k)}<0);

            nodes = bndr.Points(bndr.ConnectivityList(-f,:),:);
            n = cross(nodes(2,:)-nodes(1,:),nodes(3,:)-nodes(1,:));
            n = n/norm(n,2);
            s0= cells(i);
            e = unique(E(bis,:));
            I = 1:numel(e);
            s0Id = e==s0;
            I = [I(s0Id), I(~s0Id)];
            s = e(~s0Id);
            s1 = s(1); s2 = s(2);
            A = [(pts(s1,:)-pts(s0,:));(pts(s2,:)-pts(s0,:));n(1,:)];
            B = [V(tri(j,k),:)-pts(s0,:),(pts(s1,:)-V(tri(j,k),:)),zeros(1,3);...
                 V(tri(j,k),:)-pts(s0,:), zeros(1,3), pts(s2,:)-V(tri(j,k),:);...
                 zeros(1,9)];
            [e,I] = sort(e);
            gradC = A\B;
            gradC = [zeros(3,3*(e(1)-1))     ,gradC(:,3*I(1)-2:3*I(1)),...
                     zeros(3,3*(e(2)-e(1)-1)),gradC(:,3*I(2)-2:3*I(2)),...
                     zeros(3,3*(e(3)-e(2)-1)),gradC(:,3*I(3)-2:3*I(3)),...
                     zeros(3,numel(pts)-3*e(3))];
          case 3
            bis = symV{tri(j,k)};
            s0 = cells(i);

            e = unique(E(bis,:));
            I = 1:numel(e);
            s0Id = e==s0;
            I = [I(s0Id), I(~s0Id)];
            s = e(~s0Id);
            s1 = s(1); s2 = s(2); s3 = s(3);
            A = [(pts(s1,:)-pts(s0,:));pts(s2,:)-pts(s0,:);pts(s3,:)-pts(s0,:)];
            B = [V(tri(j,k),:)-pts(s0,:),(pts(s1,:)-V(tri(j,k),:)),zeros(1,6);...
                 V(tri(j,k),:)-pts(s0,:), zeros(1,3), pts(s2,:)-V(tri(j,k),:),zeros(1,3);...
                 V(tri(j,k),:)-pts(s0,:), zeros(1,6), pts(s3,:)-V(tri(j,k),:)];

            gradC = A\B;
            gradC = [zeros(3,3*(e(1)-1))     ,gradC(:,3*I(1)-2:3*I(1)),...
                     zeros(3,3*(e(2)-e(1)-1)),gradC(:,3*I(2)-2:3*I(2)),...
                     zeros(3,3*(e(3)-e(2)-1)),gradC(:,3*I(3)-2:3*I(3)),...
                     zeros(3,3*(e(4)-e(3)-1)),gradC(:,3*I(4)-2:3*I(4)),...
                     zeros(3,numel(pts)-3*e(4))];
          otherwise
            warning('this should not happen!')
        end
        
        g = g + dvdc*gradC;
      end
      g = g + [zeros(1,3*(cells(i)-1)),dvds,zeros(1,numel(pts)-3*cells(i))];
    end
  end

end


function [n,x0] = normalOfThreePts(p)
  assert(size(p,2)==9);
  p0 = p(:,1:3); p1 = p(:,4:6); p2 = p(:,7:9);

  n = cross(bsxfun(@minus, p1, p0),bsxfun(@minus, p2,p0));
  n = bsxfun(@rdivide, n,sqrt(sum(n.^2,2)));
  x0 = [mean(p(:,[1,4,7]),2),mean(p(:,[2,5,8]),2),mean(p(:,[3,6,9]),2)];

end



function n = matrixNormals(bndr,f,V)
  edg = isOnEdge(bndr,f,V);
  nodes = bndr.Points(bndr.ConnectivityList(f,:),:);
  vec = nodes(edg,:)-edg(circshift(edg,[1,0]),:);
  vec = vec/norm(vec,2);
  n = cross(nodes(2,:)-nodes(1,:),nodes(3,:)-nodes(1,:));
  n = n/norm(n,2);
  n(2,:) = cross(n,vec);
end


























