function [G] =restrictedVoronoiDiagram(p, bound)

dt = delaunayTriangulation(p);
E  = edges(dt);

surf = bound.freeBoundary;

edge2Vert = [surf(:,1:2),sort(surf(:,2:3),2),surf(:,[1,3])]';
edge2Vert = reshape(edge2Vert, 2,[])';
[~, IC] = ismember(edge2Vert,bound.edges,'rows');
face2edge = reshape(IC',3, [])';


s = dsearchn(dt.Points,sum(bound.Points(surf(1,:),:)/3,1));

Q = [1, s];
V = [];
symV = cell(0);
C = cell(size(p,1),1);
CT = cell(numel(C),1);
CT{s} = 1;
disp(E)
while ~isempty(Q)
    t =  Q(end,1); s = Q(end,2);
    Q = Q(1:end-1,:);
    NC = [E(:,2)==s, E(:,1)==s];
    bisect = find(any(NC,2));

    NT = findNeighbours(surf, t);    %sum(ismember(face2vert,face2vert(t,:)),2)==2;


    n = bsxfun(@minus, dt.Points(E(NC),:), dt.Points(s,:));
    n = bsxfun(@rdivide, n,sqrt(sum(n.^2,2)));
    x0 = bsxfun(@plus, dt.Points(E(NC),:), dt.Points(s,:))/2;


    symT = {-find(any(surf==surf(t,1),2)); ...
            -find(any(surf==surf(t,2),2)); ...
            -find(any(surf==surf(t,3),2))};

    
    [newVertex, symT] = clipPolygon(bound.Points(surf(t,:),:), n,x0,symT,symV,bisect);
    symV = [symV; symT];
    C{s} = [C{s}, size(V,1)+1:size(V,1)+size(newVertex,1)];
    V = [V;newVertex];
    [Q,CT] = updateQue(Q, symT, CT, E, NC, s, t);    

end

% Remove unwanted vertices (e.g. duplicates) 
[V, C] = cleanUpGrid(V, C, symV);

%% Add inner vertices, i.e. intersection of three bisections




[VInt, CInt] = dt.voronoiDiagram();
% Remove points outside domain
keep = [false;~isnan(bound.pointLocation(VInt(2:end,:)))];
% keep = [1+0.5*sin(pi*VInt(:,1))>VInt(:,3)&...
%               -1+0.5*sin(pi*VInt(:,1))<VInt(:,3)&...
%                -1<VInt(:,1)&VInt(:,1)<1& ...
%                -1<VInt(:,2)&VInt(:,2)<1;];
VInt = VInt(keep,:);
keepNum = find(keep);
newIdx = zeros(size(VInt,1),1);
newIdx(keepNum) = (1:size(VInt,1))';
CInt = cellfun(@(c) newIdx(intersect(c,keepNum))', CInt,'uniformOutput',false);
ni = size(VInt,1);
V = [VInt;V];
C = cellfun(@(cint,cext) [cint,cext+ni], CInt, C,'uniformOutput',false) ;

% for i = 1:numel(C)
%     if numel(C{i})>2
%     hull = convhulln(V(C{i},:));
%         patch('Vertices', V(C{i},:), 'faces', hull,'facecolor','y')
%     end
% end

G = voronoi2mrst(V,C, false(numel(C),1),'pebi')


end

function [V, C] = cleanUpGrid(V, C,symV)
    % Remove duplicate vertexes
    symV = cellfun(@sort, symV,'UniformOutput', false);
    symV = cellfun(@(c) num2str(c'), symV,'UniformOutput', false);
    [~,IA,IC] = unique(symV);
    VO = V;
    CO = C;
    C = cellfun(@(c) unique(IC(c)'), C,'UniformOutput',false);
    V = V(IA,:);
    
end

function [symV] = updateSym(localSym, NC, NT)
    if localSym<0
        symV = -NT(-localSym);
    else
        symV = NC(localSym);
    end
end


function [Q, CT] = updateQue(Q, symV, CT, E, NC, s, t)
    % Find possible new cells
    symV = cell2mat(symV);
    bNew = unique(symV(symV>0));
    tNew = -unique(symV(symV<0));
    for i = 1:numel(bNew)
       if isempty(CT{E(bNew(i),NC(bNew(i),:))}) || ~any(CT{E(bNew(i),NC(bNew(i),:))}==t) %New cell facet pair
           Q = [Q; t, E(bNew(i),NC(bNew(i),:))];
           CT{E(bNew(i),NC(bNew(i),:))} = [CT{E(bNew(i),NC(bNew(i),:))}, t];
       end
    end
    for i = 1:numel(tNew)
        if ~any(CT{s}==tNew(i))
           Q = [Q; tNew(i), s];
           CT{s} = [CT{s}, tNew(i)];
        end
    end
end


function NT = findNeighbours(V, t)
    VT = V(t,:);
    V  = [V(1:t-1,:);nan,nan,nan;V(t+1:end,:)];
    NT = [t;...
          find(sum(ismember(V,VT([1,2])),2)==2);...
          find(sum(ismember(V,VT([2,3])),2)==2);...
          find(sum(ismember(V,VT([3,1])),2)==2)];
end
























