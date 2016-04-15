function [p, symP] = clipPolygon(p, n, x0, symP, symV, bisector)

assert(size(n,1)==size(x0,1),'unconsistent size of n and x0');
assert(size(p,1)>2, 'A polygon needs more than 2 vertexes');


%IC = zeros(size(p,1),1);
% symP = [-4,-2,-1;
%         -3,-2,-1;
%         -4,-3,-1];
for i = 1:size(n,1)
    
    % find distance from point to plane
    d = bsxfun(@ minus, p, x0(i,:))*n(i,:)';
    
    if all(d>0)
        p  = [];
        %IC = [];
        symP = [];
        disp('wtf???')
        return
    end
    p = [p(end,:);p];
    d = [d(end);  d];
    symP = [symP{end}; symP];   

    c1 = d(1:end-1)<=0 & d(2:end)<=0;
    c2 = sign(d(1:end-1))~=sign(d(2:end));
    c3 = (c1) | (c2 & d(2:end)<0);
    

    c2 = find(c2);
    c3 = find(c3);
    
    symP = [intersectUnion(symP,c2, bisector(i)); symP(c3+1)];


    alpha = abs(d(c2))./(abs(d(c2))+abs(d(c2+1)));
    p  = [bsxfun(@times,p(c2+1,:)-p(c2,:),alpha)+p(c2,:);p(c3+1,:)];
    %IC = [i*ones(size(c2,1),1);IC(c3)];
    
    c = [c2;c3];
    [~,I] = sort(c);
    p   = p(I,:);
    %IC = IC(I);
    symP = symP(I);
end


 
end


function [B] = intersectUnion(A, idx, b)
    B = cell(size(idx,1),1);
    for i = 1:size(idx,1)
       B{i} = [intersect(A{idx(i)},A{idx(i)+1});b]; 
    end
    
end


function c = colonVec(a,b)
    c = arrayfun(@colon, a,b,'uniformOutput',false);
    c = cell2mat(c');
end