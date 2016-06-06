function [pts, removed] = removeConflictPoints3(pts,gs, pri)

[pri,I] = sort(pri);
pts = pts(I,:);
gs = gs(I);

[priNum,~, IC] = unique(pri);

removed = false(size(pts,1),1);
for i = 1:numel(priNum)-1
  low = IC==i;
  high = IC>i;
  p1 = pts(low,:);
  p2 = pts(high,:);
  g2 = gs(high);
  
  d = pdist2(p1,p2);
  rem =  any(bsxfun(@le, d,g2'),2);
  removed(rem) = true;
end

pts = pts(~removed,:);

end