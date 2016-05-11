function [G] = mirrorPebi2D(pts, bnd)
  bnd = [bnd;bnd(1,:)];
  
  t = bnd(2:end,:)- bnd(1:end-1,:);
  t = bsxfun(@rdivide, t,sqrt(sum(t.^2,2)));
  n = [t(:,2),- t(:,1)];
  
  x0 = bnd(1:end-1,:);

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