function [tout]=modify_tp(t,p,year1);

t1=t(year1);
p1=p(year1);

%linear model of the two using detrended data
t1=t1-mean(t1);
p1=p1-mean(p1);
r=polyfit(detrend(p1),detrend(t1),1);

% this is a crude estimate of the maximum amount of influence that
% precipitation variability has had on TASMAX variability
l=polyval(r,p-mean(p(year1)));

tout=t-l;

