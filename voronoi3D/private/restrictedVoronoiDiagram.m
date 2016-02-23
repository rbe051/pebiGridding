function [G] =restrictedVoronoiDiagram(p, bound)

dt = delaunayTriangulation(p);
E  = edges(dt);

face2vert = bound.freeBoundary;

edge2Vert = [face2vert(:,1:2),sort(face2vert(:,2:3),2),face2vert(:,[1,3])]';
edge2Vert = reshape(edge2Vert, 2,[])';
[~, IC] = ismember(edge2Vert,bound.edges,'rows');
face2edge = reshape(IC',3, [])';


s = dsearchn(dt.Points,sum(bound.Points(face2vert(1,:),:)/3,1));

aa = patch('Vertices', bound.Points, 'faces', face2vert,'facealpha',0.2);
hold on
str = 'e';
% for i = 1:size(p,1)
%    plot3(p(i,1),p(i,2),p(i,3),'.','markersize',30)
%    str = [str;num2str(i)];
% end
%legend(str)
Q = [1, s];
V = [];
vId = [];
C = cell(size(p,1),1);
CT = cell(numel(C),1);
CT{s} = 1;
while ~isempty(Q)
    t =  Q(end,1); s = Q(end,2);
    Q = Q(1:end-1,:);
    NC = [E(E(:,1)==s,2); E(E(:,2)==s,1)];
    
    NT = findNeighbours(face2vert, t);%sum(ismember(face2vert,face2vert(t,:)),2)==2;


    n = bsxfun(@minus, dt.Points(NC,:), dt.Points(s,:));
    n = bsxfun(@rdivide, n,sqrt(sum(n.^2,2)));
    x0 = bsxfun(@plus, dt.Points(NC,:), dt.Points(s,:))/2;
        
    a = patch('Vertices', bound.Points, 'faces', face2vert(t,:),'facealpha',0.2);
    hold on
    [newVertex, symV] = clipPolygon(bound.Points(face2vert(t,:),:), n,x0);
    
    C{s} = [C{s}, size(V,1)+1:size(V,1)+size(newVertex,1)];
    V = [V;newVertex];
    
    [Q,CT] = updateQue(Q, symV, CT, NT, NC, s, t);

    
    delete(a);
    
end

G = []

end


function [Q, CT] = updateQue(Q, symV, CT,NT, NC, s, t)
    % Find possible new cells
    sNew = unique(symV(symV>0));
    tNew = -unique(symV(symV<0));
    for i = 1:numel(sNew)
       if isempty(CT{NC(sNew(i))}) || ~any(CT{NC(sNew(i))}==t) %New cell facet pair
           Q = [Q; t, NC(sNew(i))];
           CT{NC(sNew(i))} = [CT{NC(sNew(i))}, t];
       end
    end
    for i = 1:numel(tNew)
        if ~any(CT{s}==NT(tNew(i)))
           Q = [Q; NT(tNew(i)), s];
           CT{s} = [CT{s}, NT(tNew(i))];
        end
    end
end


function NT = findNeighbours(V, t)
    VT = V(t,:);
    V  = [V(1:t-1,:);nan,nan,nan;V(t+1:end,:)];
    NT = [find(sum(ismember(V,VT([1,2])),2)==2);...
          find(sum(ismember(V,VT([2,3])),2)==2);...
          find(sum(ismember(V,VT([3,1])),2)==2)];
end
























