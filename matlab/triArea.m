% calculate area of triangles from STL
% stl vertices in 'vertices'

function [segAreas, varargout] = triArea(vertices, segments)
    % all inputs are required; 
    % vertices - an N-by-3-by-3 matrix where each 
	% row is a triangle, each column a point with (x,y,z) triplet
	% triangle 1 has vertices: A(1,1,1),(1,1,2),(1,1,3); B(1,2,1), etc
	% (triNumber,pointNum,x/y/z)
	% segments - an N-by-4 matrix with (segNumber,x,y,z); N = all 
	% skelpoints (that contribute to a segment)
	% endpoints - not needed? 
	
	% this function calls triInfo (sub-function), closestPoint.m (

areas=zeros(length(vertices));
centroids=zeros(length(vertices),3);

for i=1:length(vertices(:,1,1)); % for all STL vertices
	
	currentVertex = vertices(i,:,:); % get current 3 vertices
	[area, centroid] = triInfo(currentVertex); % call function triInfo
	
	areas(i) = area; % log current area & centroid
	centroids(i,:) = centroid;
	
	

end

