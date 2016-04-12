% point distance

function [dist] = point_dist(point1,point2)

x1=point1(1); y1=point1(2); z1=point1(3);
x2=point2(1); y2=point2(2); z2=point2(3);

dist = sqrt( (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
