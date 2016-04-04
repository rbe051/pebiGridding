function [Pts, removed] = enforceSufficientFaultCondition(Pts, CC, CR)
  I = conflictCircles(Pts, CC, CR);
  removed = false(size(Pts,1),1);
  removed(vertcat(I{:})) = true;
  Pts = Pts(~removed,:);
end
