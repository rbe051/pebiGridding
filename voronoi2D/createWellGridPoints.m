function [wellPts, wGs] = createWellGridPoints(wellLines, wellGridSize, varargin) 
% Places well grid points along wells.
%
% SYNOPSIS:
%   [wellPts, wGs] = createWellGridPoints(F,faultGridSize)
%   [...] = createWellGridPoints(..., 'Name1', Value1, 'Name2', Value2,...)
%
% Parameters:
%   wellLines       A cell of arrays. Each nx2 array in the cell contains 
%                   the piecewise linear approximation of a well. The 
%                   values must be sorted along the line, e.g., a well 
%                   consisting of two lines would be [x1,y1; x2,y2; x3,y3].
%                      .----------.-------------.
%                   (x1,y1)    (x2,y2)       (x3,y3)
%
%                   If an array has length 1 it is considered a point well.
%
%   wellGridSize    Desired distance between well points along a well.
%                   The true distance is set such that it is a factor of 
%                   the total fault length, and therefore might be slightly
%                   smaller than the desired distance.
%   
%   wfCut           - OPTIONAL.
%                   Default value array of zeros. Array of length equal the
%                   number of wells. The value equals the output of the 
%                   function [~,~, fwCut] = splitWells. The value of 
%                   element i tells if the start or end point of well i 
%                   should be removed. If the value is 1 the end point is
%                   removed. If the value is 2 the start point is removed.
%                   If the value is 3 both the starts and end point is 
%                   removed.
%
% RETURNS:
%   wellPts         Array of all generated points.
%   wGs             Distance between consecutive well points.
%
% EXAMPLE
%   wl = {[0.2,0.2;0.8,0.8]};
%   gS = 0.1;
%   wPts = createWellGridPoints(fl,gS);
%   figure(); hold on
%   plot(wl{1}(:,1), wl{1}(:,2))
%   plot(wPts(:,1),wPts(:,2),'.r','markersize',20)
%
% SEE ALSO:
%   pebiGrid, compositePebiGrid, createFaultGridPoints, splitWells, pebi.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% load options
opt   = struct('wfCut',zeros(numel(wellLines),1));
opt   = merge_options(opt,varargin{:});
wfCut = opt.wfCut;

wGs     = [];
wellPts = [];
for i = 1:numel(wellLines)  % create well points
  wellLine       = wellLines{i};
  
  if (size(wellLine,1) == 1)
      p = wellLine;
      wellSpace = wellGridSize*(1-10^-6);
  else
      [p, ~] = eqInterpret(wellLine, wellGridSize, [0;0]);
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
