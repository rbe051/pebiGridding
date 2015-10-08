function [dists] = getDistance(distance, i, j)
    %Returns the distance between point i and j
    m = 0.5*(sqrt(8*length(distance) + 1) + 1);
    assert((0<i && i<=m) && (0<j &&j<=m))
    
    if j== i
        dists = 0
    elseif i < j
        dists = distance((i-1)*(m - 0.5*i) + j - i);
    else
        dists = distance((j-1)*(m - 0.5*j) + i - j);
    end
end