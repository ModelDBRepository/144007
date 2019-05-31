function [min_dist, min_i] =findClosestPt(ref_point, pts)
% function [min_dist, min_i]=findClosestPt(ref_point, pts)
%
% This function finds the closest point in the n-by-2 vector of pts to the
% ref_point.  Returns the distance and the index of the closest point.
% The method is very simple, and maybe not the fastest for matlab, but 
% when I tried it by making pair-wise distance matrices it seemed to take forever.

min_dist = -1;
min_i = [];
for i=1:size(pts,1)
    dist = sqrt(sum((pts(i,:)-ref_point).^2));
    if (min_dist == -1 && dist ~= 0)
        min_dist = dist;
        min_i = i;
    elseif (dist ~= 0)
        if (dist < min_dist)
            min_dist = dist;
            min_i = i;
        end  
    end
end
