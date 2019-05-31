function [X, Ymean, Ysd] = resampleVectorEvenly(x, y, sampling, nbins)
% function [X, Ymean, Ysd] = resampleVectorEvenly(x, y, sampling, nbins)

x=x(:);
y=y(:);
if strcmp(sampling, 'linear')
    bin_edges = linspace(min(x), max(x), nbins);
elseif strcmp(sampling, 'log')
    bin_edges = logspace(log10(min(x)), log10(max(x)), nbins);
else
    error('sampling needs to be either LINEAR or LOG (lowercase)');
end
X = NaN*zeros(length(bin_edges)-1,1);
Ymean = NaN*zeros(length(bin_edges)-1,1);
Ysd = NaN*zeros(length(bin_edges)-1,1);
for i=1:(length(bin_edges)-1)
    X(i) = bin_edges(i) + (bin_edges(i+1)-bin_edges(i))/2;
    si = (x >= bin_edges(i) & x < bin_edges(i+1));
    if sum(si)
        Ymean(i) = mean(y(si));
        Ysd(i) = std(y(si));
    end
end
 
X=X(:);
Ymean=Ymean(:);
Ysd = Ysd(:);