function [P1, removed] = removeConflictPoints2(P1,P2,TOL)
  % Remove conflict points
  % Arguments:
  %   P1              2*n array of points
  %   P2              2*m array of points
  %   TOL             array of length n containing the minimum allowed
  %                   distance around P2
  %
  % Return:
  %   P1              array of length <=n of points that were not removed
  %   removed         logical array of length n that is true for points
  %                   that were removed
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  TOL = repmat(TOL',size(P1,1),1);
  removed = any(pdist2(P1,P2)<TOL,2);
  P1      = P1(~removed,:);

end