function wo = weibull_func3(params, x)
% function wo = weibull_func(params, x)
%
% This is a function for the fitting of weibull CDF parameters to data.  Suitable
% to use for lsqcurvefit.

wo = params(1) .* wblcdf(x, params(2), params(3));

