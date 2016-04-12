% scale coordinates (assumes same orientation)
function [scaled] = scaleCoords(unscaled, reference);
% [scaled] = scaleCoords(unscaled, reference)
% input is unscale (coords TO BE scaled)
% and reference (will not be scaled)
% must be in Nx3 tuples for both

unscaled_max_x = max(unscaled(:,1));
unscaled_max_y = max(unscaled(:,2));
unscaled_max_z = max(unscaled(:,3));
unscaled_min_x = min(unscaled(:,1));
unscaled_min_y = min(unscaled(:,2));
unscaled_min_z = min(unscaled(:,3));

unscaled_diffx = unscaled_max_x - unscaled_min_x;
unscaled_diffy = unscaled_max_y - unscaled_min_y;
unscaled_diffz = unscaled_max_z - unscaled_min_z;

reference_max_x = max(reference(:,1));
reference_max_y = max(reference(:,2));
reference_max_z = max(reference(:,3));
reference_min_x = min(reference(:,1));
reference_min_y = min(reference(:,2));
reference_min_z = min(reference(:,3));

reference_diffx = reference_max_x - reference_min_x;
reference_diffy = reference_max_y - reference_min_y;
reference_diffz = reference_max_z - reference_min_z;

xfact = unscaled_diffx / reference_diffx;
yfact = unscaled_diffy / reference_diffy;
zfact = unscaled_diffz / reference_diffz;

scaled = zeros(size(unscaled));
scaled(:,1) = unscaled(:,1) / xfact;
scaled(:,2) = unscaled(:,2) / yfact;
scaled(:,3) = unscaled(:,3) / zfact;
