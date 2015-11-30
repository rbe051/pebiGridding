function [p] = createFaultLine(line, fh, lineDist, varargin)

dptol=.001; Fscale=1.2; deltat=.1; geps=.001*lineDist;

% 1. Create initial points, equally distributed.
p = eqInterpret(line, lineDist);
if size(p,1)<3, return; end

N=size(p,1);                                         % Number of points N

count=0;
clf,view(2),axis equal,axis off
while 1 & count<1000
  count=count+1;
  % 2. Calculate distances, and wanted distances
  d = sqrt(sum((p(1:end-1,:) - p(2:end,:)).^2,2));
  pmid = (p(1:end-1,:) + p(2:end,:))/2;
  dw = lineDist*fh(pmid,varargin{:});
  
  % 3. Move points based on bar desired length
  F = dw - d;                                       % Bar forces (scalars)
  d(count) = d(count) + F(count);                   % New bar lengths
  p = interpLine(line,d);                           % Update node positions
  plot(p(:,1), p(:,2), '.')
  % 4. Termination criterion: All bar length is close to wanted length(scaled)
  if count>=size(p,1)-1, break; end
end

if count == 1000
    warning('DistMesh did not converge in 1000 iterations.')
end

end