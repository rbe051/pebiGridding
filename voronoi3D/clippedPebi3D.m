function [G,dt,symV] = clippedPebi3D(p, bound)

dt = delaunayTriangulation(p);
E  = edges(dt);

dtB.ConnectivityList = bound.freeBoundary;
dtB.Points = bound.Points;


% Clip grid against boundary
[Vext,Cext,symV] = clipGrid(dt,dtB);


%% Add inner vertices, i.e. intersection of three bisections




[VInt, CInt] = voronoin(p);%dt.voronoiDiagram();
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
V = [VInt;Vext];
symV = [cell(size(VInt,1),1);symV];
C = cellfun(@(cint,cext) [cint,cext+ni], CInt, Cext,'uniformOutput',false) ;


[V,IA,IC] = uniquetol(V,100*eps,'byRows',true);
symV = symV(IA);
C = cellfun(@(c) unique(IC(c))', C,'UniformOutput',false);

% color = get(gca,'ColorOrder');
% nc = size(color,1);
% for i = 1:numel(C)
% if numel(C{i})>2
% hull = convhulln(V(C{i},:));
% patch('Vertices', V(C{i},:), 'faces', hull,'facecolor',color(mod(i,nc)+1,:),'facealpha',0.3)
% a = plot3(V(C{i},1),V(C{i},2),V(C{i},3),'.','markersize',20);
% delete(a)
% end
% end
G = voronoi2mrst(V,C, false(numel(C),1),'pebi');

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
























