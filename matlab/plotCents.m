function [surface] = plotCents (vert, type, plot)
% inputs: - vertices - expected as triangles where triangle i is given by
% vert(i,:,1) vert(i,:,2) and vert(i,:,3); 
% - if so then type = 'tri'
% - if plot = 1, plot is made
% outputs: -surface, an Mx4 matrix where each triangle M has a tuple
% (x,y,z) in (M,1:3) and an area in (M,4)

if type == 'tri';
    
    [m,n,p] = size(vert);
        
    % create matrix surface
    surface = zeros(m,4);
    
    ProgressBar('CenterPoints',m);
    
    for i = 1:m;
        
        % create new temporary vertices
        currentvert = zeros(3,3);
       % currentvert = vert(i,:,:);
        point1 = vert(i,:,1);
        point2 = vert(i,:,2);
        point3 = vert(i,:,3);
        currentvert(:,1) = point1;
        currentvert(:,2) = point2;
        currentvert(:,3) = point3;
        
        % get centroid
        [centroid, area] = triInfo(currentvert);
%         centroids(i,:) = centroid;
%         areas(i) = area;
        surface(i,1:3) = centroid;
        surface(i,4) = area;
        
        ProgressBar('CenterPoints');
        
    end
    
    
else
end



if plot == 1;
    fprintf('Plotting scatter of surface points...')
    scatter3(surface(:,1),surface(:,2),surface(:,3),'k.');
else
end

    