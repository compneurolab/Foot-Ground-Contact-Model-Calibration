function [fitdata,fitdatap] = SplineSmooth(time,data,param)

% Smooth input curve and output smoothed curve and its first derivative

% Determine number of rows of data
mrows = length(time);

% Set up curve-fitting problem
ftype = fittype('smoothingspline' );
fopts = fitoptions('method','SmoothingSpline','SmoothingParam',param);

% Perform spline smoothing
fitmodel = fit(time,data,ftype,fopts);
pp = fitmodel.p;
fitdata = ppval(pp,time);

% Calculate first derivative of smoothed curve
pp_p = pp;
pp_p.coefs(:,1) = 0;
pp_p.coefs(:,2) = 3*pp.coefs(:,1);
pp_p.coefs(:,3) = 2*pp.coefs(:,2);
pp_p.coefs(:,4) = 1*pp.coefs(:,3);
fitdatap = ppval(pp_p,time);
