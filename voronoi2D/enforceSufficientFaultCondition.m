function [Pts, removed] = enforceSufficientFaultCondition(Pts, CC, CR)
% Remove any points that innside given circles
%
% SYNOPSIS:
%   [Pts, removed] = enforceSufficientFaultCondition(Pts, CC, CR)
%
% PARAMETERS;
%   Pts             A nx2 array of possible conflict points. If a point is
%                   inside any of the circles given by CC and CR it will be
%                   removed.
%   CC              A mx2 array of the center coordinates of the circles 
%
%   CR              A mx1 array of the circle radii
%
% RETURNS:
%   Pts             All points outside the circles
%
%   removed         A nx1 logical array that is true for any Pts innside a
%                   circle
%
% EXAMPLE
%   [X,Y] = meshgrid(1:10,1:10);
%   pts = [X(:),Y(:)];
%   cc = [5,5;2,2];
%   cr = [3,2];
%   ptsRem = enforceSufficientFaultCondition(pts,cc,cr);
%   figure(); hold on
%   plot(pts(:,1),pts(:,2),'o')
%   plot(ptsRem(:,1),ptsRem(:,2),'.')
%   theta = linspace(0,2*pi)'
%   for i = 1:size(cc)
%   X = repmat(cc(i,:),100,1)+repmat(cr(i),100,2).*[cos(theta),sin(theta)]
%   plot(X(:,1), X(:,2))
%   end
%
% SEE ALSO:
%   splitWells, createFaultGridPoints, compositePebiGrid, pebiGrid

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
  I = conflictCircles(Pts, CC, CR);
  removed = false(size(Pts,1),1);
  removed(vertcat(I{:})) = true;
  Pts = Pts(~removed,:);
end
