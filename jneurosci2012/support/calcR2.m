function rsq = calcR2(ydata, yfit) 

yresid = ydata - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(ydata)-1) * var(ydata);
rsq = 1 - SSresid/SStotal;

