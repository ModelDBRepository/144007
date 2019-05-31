%
% For plotting a nice overlay with SDs etc. - pre AND post data should be given
% if possible
%
function plot_error (time, pre_mean, pre_sd, post_mean, post_sd, pre_color, pre_patch, post_color, post_patch)
  % First plot the SDs
  [x_err_poly, y_err_poly] = get_sem_poly(time, pre_mean, pre_sd);
  patch(x_err_poly,y_err_poly, pre_patch, 'EdgeColor', 'none');
  if (length(post_mean) > 0)
    [x_err_poly, y_err_poly] = get_sem_poly(time, post_mean, post_sd);
    patch(x_err_poly,y_err_poly, post_patch, 'EdgeColor', 'none');
  end
  
  % Then the lines around SDs
%  plot(time, pre_mean + pre_sd, [pre_color '-']);
%  plot(time, pre_mean - pre_sd, [pre_color '-']);
%  plot(time, post_mean + post_sd, [post_color '-']);
%  plot(time, post_mean - post_sd, [post_color '-']);
  
  % And finally the main line
  plot(time, pre_mean, 'Color', pre_color, 'Linewidth', 2);
  if (length(post_mean) > 0)
    plot(time, post_mean, 'Color', post_color, 'Linewidth', 2);
  end

  
  
%
% Returns your data as a polygon for bounding y at each x value with +/- y_off
%  for real nice SEM/SD plots
%
function [ret_x, ret_y] = get_sem_poly(x, y, y_off);
  ret_x = zeros(2*length(x),1);
  ret_y = zeros(2*length(x),1);
  
  l = length(x);
  
  for i=1:length(ret_x)
    if (i < l)
      ret_x(i) = x(i);
      ret_y(i) = y(i) + y_off(i);
    elseif (i == l || i == l+1)
      ret_x(i) = x(l)+2;
      ret_y(i) = y(l) + y_off(l);
      if ( i == l + 1) ; ret_y(i) = y(l) - y_off(l); end
    else
      ret_x(i) = x(2*l - i + 1);
      ret_y(i) = y(2*l - i + 1) - y_off(2*l - i + 1);
    end
  end
