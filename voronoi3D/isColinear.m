function [id] = isColinear(pts)
    id = rank(bsxfun(@minus,pts(2:end,:),pts(1,:)),1e-10)<2;
end