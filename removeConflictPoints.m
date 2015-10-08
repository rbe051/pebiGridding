function [Pts, removed] = removeConflictPoints(Pts, gridSpacing, priIndex)
    Ic = 1:size(Pts, 1);
    ptsToClose = Pts;
    removed = zeros(size(Pts, 1), 1);
    
    distance = pdist(ptsToClose)';
    dlt = distLessThan(distance, gridSpacing(Ic));
    Ic = findToClose(dlt);
    
    while length(Ic)>1
        sumToClose = sumToClosePts(dlt);
        sumToClose = sumToClose(find(sumToClose));
        [~, Is] = sort(sumToClose,'descend');
        [~, Ii ] = sort(priIndex(Ic(Is)), 'ascend');

        removePoint = Ic(Is(Ii(1)));
        removed(removePoint) = 1;        
        Ic = Ic(Ic~=removePoint);
        ptsToClose = Pts(Ic,:);
        
        if ptsToClose ==1
            continue
        end
        distance = pdist(ptsToClose)';
        dlt = distLessThan(distance, gridSpacing(Ic));
        Ic = Ic(findToClose(dlt));
    end
    Pts = Pts(~removed,:);
end


function [arr] = distLessThan(distance, b)
    n = length(distance);
    [i,j] = arrToMat(1:n, n);
    arr = distance < max(b(i), b(j));
end

function [pts] = sumToClosePts(arr)
    n = length(arr);
    m = ceil(sqrt(2*n)); % = 0.5 + 0.5sqrt(1+8n) for n > 0
    pts = zeros(m,1);
    for i = 1:m
        k1 = matToArr(i+1:m, i, m);
        k2 = matToArr(i,1:i-1, m);
        pts(i) = sum(arr(k1)) + sum(arr(k2)); 
    end
end

function indexes = findToClose(arr)
    n = length(arr);
    k = find(arr);
    [i, j] = arrToMat(k, n);
    indexes = unique([i ; j]);
end

function [k] = matToArr(i,j, m)
    assert(all(abs(j)) && all(abs(i-j)) && all(abs(1+m-i)));
    k = 1 + (j-1)*m - (j-1).*j/2 + i-j - 1;
end

function [i, j] = arrToMat(k, n)
    m = ceil(sqrt(2*n));
    j = ceil((2*m-1)/2 - 0.5*sqrt((2*m-1)^2 - 8*k));
    i = k + j - (j-1).*(m-j/2);
end