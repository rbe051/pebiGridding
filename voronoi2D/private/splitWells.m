function [wellLines, wCut, wfCut] = splitWells(faultLines, wellLines)
tempWell = cell(0);
tmpCut    = cell(0);

for i = 1:numel(wellLines)
  [tempWell{i},~,tmpCut{i}]  = splitLines(wellLines(i),wellLines([1:i-1,i+1:end]));
end
tempWell = horzcat(tempWell{:});
tmpCut      = vertcat(tmpCut{:});
[wellLines, cutId, wfCut] = splitLines(tempWell, faultLines);

wCut = [];
for i = 1:size(tmpCut,1)
  switch tmpCut(i)
    case 0
      wCut = [wCut; zeros(sum(cutId(i,:)),1); 0];
    case 1
      wCut = [wCut; zeros(sum(cutId(i,:)),1); 1];
    case 2
      wCut = [wCut; 2; zeros(sum(cutId(i,:)),1)];
    case 3 
      if sum(cutId(i,:))
        wCut = [wCut;2; zeros(sum(cutId(i,:)-1),1);1];
      else
        wCut =  [wCut;3];
      end
  end
end

end
