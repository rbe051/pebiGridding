function [faultLines, fCut, fwCut] = splitFaults(faultLines, wellLines)

tempFault = cell(0);
tmpCut    = cell(0);

for i = 1:numel(faultLines)
  [tempFault{i},~,tmpCut{i}] = splitLines(faultLines(i),faultLines([1:i-1,i+1:end]));
end

tempFault = horzcat(tempFault{:});
tmpCut      = vertcat(tmpCut{:});
[faultLines, cutId,fwCut] = splitLines(tempFault, wellLines);

fCut = [];
for i = 1:size(tmpCut,1)
  switch tmpCut(i)
    case 0
      fCut = [fCut; zeros(sum(cutId(i,:)),1); 0];
    case 1
      fCut = [fCut; zeros(sum(cutId(i,:)),1); 1];
    case 2
      fCut = [fCut; 2; zeros(sum(cutId(i,:)),1)];
    case 3 
      if sum(cutId(i,:))
        fCut = [fCut;2; zeros(sum(cutId(i,:)-1),1);1];
      else
        fCut =  [fCut;3];
      end
  end
end

end