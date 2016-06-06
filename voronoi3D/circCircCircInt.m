function [p,offSet] = circCircCircInt(CC1, CR1,CC2,CR2,CC3,CR3)

if isempty(CC1) || isempty(CC2)
  p = [];
  return
end
assert(all(size(CC1)==size(CC2) & size(CC1)==size(CC3)),...
  'Circle centers must have the same size');
assert(size(CC1,1)==size(CR1,1) ...
    && size(CC2,1)==size(CR2,1) ...
    && size(CC3,1)==size(CR3,1),...
    'The number of circle centers must equal the number of radiuses');

dim = size(CC1,2);


p1 = CC2 - CC1;
dx = sqrt(sum((p1).^2,2));               % Distance between c1 & c2f
X = (dx.^2 - CR2.^2 + CR1.^2)./(2*dx);          % Mid-Point
%faultOffsetx = sqrt(CR1.^2 - bisectPntx.^2);    % Pythagoras
nx = bsxfun(@rdivide,p1,dx);               % Unit vector

%dy = sqrt(sum((CC1 - CC3).^2,2));               % Distance between c1 & c2f
p2 = CC3 - CC1;
midX = dot(nx, p2,2);
ny   = p2 - bsxfun(@times,midX,nx);
ny   = bsxfun(@rdivide, ny, sqrt(sum(ny.^2,2)));
midY = dot(ny,p2,2);
Y = (CR1.^2 - CR3.^2 + midY.^2 + midX.^2)./(2*midY) - midX./midY.*X;

nz = cross(nx,ny);
Z  = sqrt(CR1.^2 - X.^2-Y.^2);        % Pythagoras

% Set right left and right intersection points
top = CC1 + bsxfun(@times, X, nx)  ...
          + bsxfun(@times, Y, ny)  ...
          + bsxfun(@times, Z, nz);
bot = top -2*bsxfun(@times, Z, nz);

% Put result together
p = reshape([top,bot]',dim,[])';
offSet = reshape([Z,Z]',1,[])';
% figure()
% hold on
% [X,Y,Z] = sphere;
% for i = 1:size(CC1,1)
%   a = plot3(p(2*i-1:2*i,1),p(2*i-1:2*i,2), p(2*i-1:2*i,3),'.','markersize',20)
%   b = surf(CC1(i,1) + X*CR1(i),CC1(i,2) + Y*CR1(i),CC1(i,3) + Z*CR1(i))
%   c = surf(CC2(i,1) + X*CR2(i),CC2(i,2) + Y*CR2(i),CC2(i,3) + Z*CR2(i))
%   d =surf(CC3(i,1) + X*CR1(i),CC3(i,2) + Y*CR3(i),CC1(i,3) + Z*CR3(i))
%   delete(a)
%   delete(b)
%   delete(c)
%   delete(d)
% end

end