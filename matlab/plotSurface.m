function plotSurface(points)
% usage: plotSurface(points)
% Uses meshgrid, griddata and surf to plot scattered points as a surface

xs = points(:,1);
ys = points(:,2);
zs = points(:,3);

dx = min(diff(xs)) * 0.5;
dy = min(diff(ys)) * 0.5;
dz = min(diff(zs)) * 0.5;

x_edge = [min(xs) - dx : dx : max(xs) + dx];
lenx = length(x_edge);
y_edge = [min(ys) - dy : dy : max(ys) + dy];
leny = length(y_edge);
z_edge = [min(zs) - dz : dz : max(zs) + dz];

[Xs, Ys] = meshgrid(x_edge, y_edge);

Zs = griddata(xs,ys,zs,Xs,Ys);

surf(Xs,Ys, Zs)