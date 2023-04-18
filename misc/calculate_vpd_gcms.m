function dataout=calculate_vpd_gcms(rhmean,tmin,tmax)



%VPD:
%Calculate Vapor Pressure Deficit (in KPa, Kilopaschals) in R according to equations of Allen et al. (1998) for Example 5 using mean monthly minimum (tmin) and maximum (tmax) temperatures (Celsius) and mean monthly pecent relative humidity (rh, 100 = 100%): [http://www.fao.org/docrep/x0490e/x0490e07.htm#chapter 3 meteorological data]
%Declare in R values of three required environmental variables from Example 5 of Allen et al. (1998)

esmn=0.6108 * exp((17.27 * tmin) ./ (tmin + 237.3));
esmx=0.6108 * exp((17.27 * tmax) ./ (tmax + 237.3));
esm= (esmn + esmx)/2;
ea= (rhmean/100) .* esm;
dataout=esm - ea;

% Calculate saturation vapor pressure for mean minimum monthly temperature (esmn)
% 
% get.esmn <- function(tmin){
%   esmn <- .6108 * exp((17.27 * tmin) / (tmin + 237.3))
%   return(esmn)
% }
% 
% Calculate saturation vapor pressure for mean maximum monthly temperature (esmx)
% 
% get.esmx <- function(tmax){
%   esmx <- .6108 * exp((17.27 * tmax) / (tmax + 237.3))
%   return(esmx)
% }
% 
% Calculate mean saturation vapor pressure (esm)
% 
% get.esm <- function(tmin, tmax){
%   esmn <- get.esmn(tmin)
%   esmx <- get.esmx(tmax)
%   esm <- (esmn + esmx)/2
%   return(esm)
% }
% 
% Calculate actual vapor pressure (ea)
% 
% get.ea <- function(rh, tmin, tmax){
%   esm <- get.esm(tmin, tmax)
%   ea <- (rh/100) * esm
%   return(ea)
% }
% 
% Calculate vapor pressure deficit (vpd = esm - ea; getting esm and ea functions)
% 
% get.vpd <- function(rh, tmin, tmax){
%   esm <- get.esm(tmin, tmax)
%   ea <- get.ea(rh, tmin, tmax)
%   vpd <- esm - ea
%   return(vpd)
% }
% 
% Check variable values and results
% 
% esmn <- get.esmn(tmin)
% esmx <- get.esmx(tmax)
% esm <- get.esm(tmin, tmax)
% ea <- get.ea(rh)
% vpd <- get.vpd(rh, tmin, tmax)
% 
% Define temp (mean temperature), tmin and tmax for plotting relationship of esm and vpd
% 
% get.temp <- function(tmin, tmax){
%   temp <- (tmin + tmax)/2
%   return(temp)
% }
% tmin <- -40:20
% tmax <- -20:40
% 
% Plot relationship of esm and vpd to mean temperature
% 
% plot(get.temp(tmin,tmax), get.esm(tmin, tmax), type = "l", xlab = "Temp (C)", ylab = "esm (black) or vpd (red) (kPa)")
% lines(temp, get.vpd(50, tmin, tmax), col = "red")