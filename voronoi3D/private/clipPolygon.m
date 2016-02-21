function [p] = clipPolygon(p, n, x0)

assert(size(n,1)==size(x0,1),'unconsistent size of n and x0');
assert(size(p,1)>2, 'A polygon needs more than 2 vertexes');


for i = 1:size(n,1)
    % find distance from point to plane
    d = bsxfun(@ minus, p,x0(i,:))*n(i,:)';
    p = calculateNewVertex(p,d);
end

end


function [newVertex] = calculateNewVertex(p, d)
    p = [p(end,:); p];
    d = [d(end)  ; d];
    c1 = d(1:end-1)>=0 & d(2:end)>=0;
    c2 = sign(d(1:end-1))~=sign(d(2:end));
    c3 = c2 & d(2:end)>0;
    
    c2 = find(c2);
        
    alpha = abs(d(c2))./(abs(d(c2))+abs(d(c2+1)));
    c1 = find(c1|c3);
    
    p = [bsxfun(@times,p(c2+1,:)-p(c2,:),alpha)+p(c2,:);p(c1+1,:)];

    c = [c2;c1];
    [~,I] = sort(c);
    newVertex = p(I,:);
end