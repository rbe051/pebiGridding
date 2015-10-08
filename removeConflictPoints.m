function [Pts, removed] = removeConflictPoints(Pts, gridSpacing, priIndex)
    distance = pdist(Pts)';
    rm = distLessThan(distance, gridSpacing);   
    
    rm = squareform(rm);
    One = ones(size(rm,1), 1);
    sumShort = rm*One;
    toClose = find(sumShort,1);
    
    [~, Is] = sort(sumShort,'descend');
    [~, Ii ] = sort(priIndex(Is), 'ascend');
    
    removed = zeros(size(Pts, 1), 1);
    
    while ~isempty(toClose)
        for i = 1:length(sumShort)
            if sum(toClose == Is(Ii(i)))
                removePoint = Is(Ii(i));
                removed(removePoint) = 1;
                rm(removePoint,:) = 0;
                rm(:, removePoint) = 0;
                sumShort = rm*One;
                toClose = find(sumShort);
                [~, Is] = sort(sumShort, 'descend');
                [~, Ii] = sort(priIndex(Is), 'ascend');
                break
            end
        end    
    end
    Pts = Pts(~removed,:);
end


function [arr] = distLessThan(distance, b)
    n = length(distance);
    m = length(b);
    k = 1:n;
    j = ceil((2*m-1)/2 - 0.5*sqrt((2*m-1)^2 - 8*k));
    i = k + j - (j-1).*(m-j/2);
    arr = distance < max(b(i), b(j));
end