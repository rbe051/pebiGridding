function [sWells, wCut, wfCut] = splitWells(faultLines, wellLines)
% Split well lines on every well-well and fault-well intersection
%
% SYNOPSIS:
%   [sFaults, fCut, fwCut] = splitWells(faultLines, wellLines)
%
% PARAMETERS;
%   faultLines      A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear approximation of a fault. The 
%                   values must be sorted along the line, e.g., a fault 
%                   consisting of two lines would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%   wellLines       A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear approximation of a well. The 
%                   values must be sorted along the line, e.g., a well 
%                   consisting of two lines would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%
%                   If an array has length 1 it is considered a point well.
% RETURNS:
%   sWells          A cell of arrays. The arrays are the cut faults and 
%                   does not contain any intersections with faultlines or 
%                   wellLines, except possible at the start and/or end 
%                   points.
%
%   wCut            Array of length equal sFaults. The value of element
%                   i tells if the returned fault i has an intersection 
%                   with faultLines. If the value is 0, the fault has no 
%                   intersections. If the value is 1 it has an intersection 
%                   at the end point. If the value is 2 it has an 
%                   intersection at the start point. If the value is 3 it 
%                   has an intersection at both the start and end point.
%
%   fwCut           Array of length equal sFaults. The value of element
%                   i tells if the returned fault i has an intersection 
%                   with wellLines. If the value is 0, the fault has no 
%                   intersections. If the value is 1 it has an intersection 
%                   at the end point. If the value is 2 it has an 
%                   intersection at the start point. If the value is 3 it 
%                   has an intersection at both the start and end point.
%
% EXAMPLE
%   fl = {[0.2,0.2;0.8,0.8], [0.2,0.5;0.8,0.5]};
%   wl = {[0.3,0.8;0.3,0.2]};
%   gS = 0.1;
%   fl = splitFaults(fl,wl);
%   figure(); hold on
%   for i = 1:numel(fl)
%     plot(fl{i}(:,1), fl{i}(:,2))
%   end
%
% SEE ALSO:
%   splitWells, createFaultGridPoints, compositePebiGrid, pebiGrid

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


tempWell  = cell(0);
tmpCut    = cell(0);

for i = 1:numel(wellLines)
  [tempWell{i},~,tmpCut{i}]  = splitLines(wellLines(i),wellLines([1:i-1,i+1:end]));
end
tempWell = horzcat(tempWell{:});
tmpCut      = vertcat(tmpCut{:});
[sWells, cutId, wfCut] = splitLines(tempWell, faultLines);

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
        wCut = [wCut;3];
      end
  end
end

end
