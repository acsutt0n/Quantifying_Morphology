
function [centroid, varargout] = triInfo(vertices)
% input is matrix vertices: (3,3) = (1/2/3,xyz)
% sideLength = sqrt( (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 )
% make sure coordinates are: (triNum,pointNum,x/y/z)
%   output: 1 - centroid (x,y,z) of triangle
%           2 - area of triangle

    side1 = sqrt( (vertices(2,1)-vertices(1,1))^2 + ...
        (vertices(2,2)-vertices(1,2))^2 + ...
        (vertices(2,3)-vertices(1,3))^2 );
    side2 = sqrt( (vertices(3,1)-vertices(2,1))^2 + ...
        (vertices(3,2)-vertices(2,2))^2 + ...
        (vertices(3,3)-vertices(2,3))^2 );
    side3 = sqrt( (vertices(1,1)-vertices(3,1))^2 + ...
        (vertices(1,2)-vertices(3,2))^2 + ...
        (vertices(1,3)-vertices(3,3))^2 );

if nargout >= 2;
%     varargout{1} = area;
    s = 0.5*(side1+side2+side3);
    varargout{1} = sqrt( s*(s-side1)*(s-side2)*(s-side3) );
else
end
    
    midx=(1/3)*( vertices(1,1)+vertices(1,2)+vertices(1,3) );
    midy=(1/3)*( vertices(2,1)+vertices(2,2)+vertices(2,3) );
    midz=(1/3)*( vertices(3,1)+vertices(3,2)+vertices(3,3) );
    centroid = [midx midy midz];
