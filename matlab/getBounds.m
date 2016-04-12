function [bounds] = getBounds(matrix);
% [bounds] = getBounds(matrix);
% get the max and min of coordinates in all dimesions
% input is a matlab matrix variable
% output is nested vars

[m,cols] = size(matrix);

bounds.xmin = min(matrix(:,1));
bounds.xmax = max(matrix(:,1));
bounds.ymin = min(matrix(:,2));
bounds.ymax = max(matrix(:,2));

if cols == 3
    bounds.zmin = min(matrix(:,3));
    bounds.zmax = max(matrix(:,3));
else
end

if cols >= 4
    error('only returning first 3 dimensions');
else
end
end
