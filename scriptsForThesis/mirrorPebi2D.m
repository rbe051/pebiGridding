function [G] = mirrorPebi2D(pts, n, x0)
  k = size(pts,1);
  mirPts = zeros(k*size(n,1),2);
  for i = 1:size(n,1)
    mirPts(i*k-k+1:i*k,:) = pts - 2*bsxfun(@minus, pts, x0(i,:))*n(i,:)'*n(i,:);
  end

  Gt = triangleGrid([pts;mirPts]);
  G = pebi(Gt);

    rem = false(G.cells.num,1);
  rem(k+1:end) = true;

  G = removeCells(G, rem);


end