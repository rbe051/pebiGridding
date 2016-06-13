function s = signDistFuncPolygon(p, varargin)
bdr = varargin{2};
s = p_poly_dist(p(:,1),p(:,2),bdr(:,1),bdr(:,2));
end