function [x] = polarToCart(theta, r)
 x = [r.*cos(theta), r.*sin(theta)];
end