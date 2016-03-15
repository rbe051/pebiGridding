function [id] = isColinear(pts)
    id = rank(bsxfun(@minus,pts(2:end,:),pts(1,:)),50*eps)<2;
end