function [lum_t, RFweights] = photRF(center, width, stim_x, stim_y, stim_movie)
% function [lum_t, RFweights] = photRF(center, width, stim)
%
% This function takes a stimulus input and and position and width,
% outputting a luminance experienced as a function of time, in the brightness 
% units given as the input.
% The center is the center of the gaussian distribution that is
% the RF.  The width is a parameter that gives the full width of the
% receptive field.  The sigma parameter of the gaussian is given by the
% SD = width/4 (since 4*sigma is the approximate 95% confidence interval 
% for a gaussian).  These parameters must be given relative to the stimulus 
% movie (ie, width and position must correspond to pixels in the movie).

% figure the RF of the photoreceptor
x_range = center(1)-width:center(1)+width;
y_range = center(2)-width:center(2)+width;
x_siz = length(x_range);
y_siz = length(y_range);
xc = floor(x_siz/2)+1;
yc = floor(y_siz/2)+1;
[X,Y] = meshgrid(x_range, y_range);
RFweights = gauss2d(center, (width/4)^2, X, Y); %mu and sigma^2 parameters

t_len = size(stim_movie,3);
lum_t = zeros(t_len, 1);
if (1==0) % look at the RF
    surf(X,Y,RFweights); colorbar;
end
% need to align that RF with the stimulus
x_offset = find(stim_x >= x_range(1),1,'first');
y_offset = find(stim_y >= y_range(1), 1, 'first');
if (~isempty(x_offset) && ~isempty(y_offset))
    xend = min(size(stim_movie,1), x_offset+x_siz) - 1; %In case the RF goes outside the stimulus, take the min
    yend = min(size(stim_movie,2), y_offset+y_siz) - 1;

    for i=1:t_len
        stim_section = stim_movie(x_offset:xend, y_offset:yend,i);
        stim_section = double(stim_section);
        weighted = RFweights(1:size(stim_section,1), 1:size(stim_section,2)) .* stim_section;
        lum_t(i) = sum(sum(weighted));
    end
end


% ------------------------------------------------------------
function weights = gauss2d(center, sigma2, pos_x, pos_y)
% function weights = gauss2d(center, sigma2, pos_x, pos_y)
% Returns the values of a 2D gaussian, where the center (X,Y) coordinates are given in 
% CENTER, and the variance of both dimensions are given by SIGMA2.  The positions for
% which to return values are specified in POS_X and POS_Y.

center_x = center(1);
center_y = center(2);
weights = 1/(2*pi* sigma2).*exp(-1*((pos_x-center_x).^2+(pos_y-center_y).^2)./(2* sigma2));