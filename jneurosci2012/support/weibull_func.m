function wo = weibull_func(params, x)
% function wo = weibull_func(params, x)
%
% This is a function for the fitting of weibull CDF parameters to data.  Suitable
% to use for lsqcurvefit.

wo =  wblcdf(x, params(1), params(2));

