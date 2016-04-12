function scatterSurface(pts, Alpha)
%
% plot 3d scatter points as a surface (ellipsoid)
% alpa = how transparent the surface is (optional) 
%

if nargin > 1
  a = Alpha;
else
  a = 0.05; % default 
end

x = pts(:,1);
y = pts(:,2);
z = pts(:,3);

% delaunay triangulation
tri = delaunay(x,y);
[r, c] = size(tri);
fprintf('%i triangles \n', r);

h = trisurf(tri, x, y, z);
axis equal;
lighting phong;
shading interp;
axis off;
alpha(a);
