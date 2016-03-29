function [wellPts, wGs] = createWellGridPoints(wellLines, wellDs, wfCut) 
    
wGs     = [];
wellPts = [];
for i = 1:numel(wellLines)  % create well points
  wellLine       = wellLines{i};
  
  if (size(wellLine,1) == 1)
      p = wellLine;
      wellSpace = wellDs*(1-10^-6);
  else
      [p, ~] = eqInterpret(wellLine, wellDs, [0;0]);
      wellSpace = sqrt(sum(diff(p,1,1).^2,2));
      wellSpace = [wellSpace(1);wellSpace;wellSpace(end)];
      wellSpace = min([wellSpace(1:end-1), wellSpace(2:end)],[],2);
  end
  
  keep = 1:size(p,1);
  switch wfCut(i)
    case 1
      keep = keep(1:end-1);
    case 2
      keep = keep(2:end);
    case 3
      keep = keep(2:end-1);
  end
  wGs     = [wGs; wellSpace(keep)];
  wellPts = [wellPts;p(keep,:)];
end

[wellPts,IA] = uniquetol(wellPts,50*eps,'byRows',true);
wGs          = wGs(IA);

end
