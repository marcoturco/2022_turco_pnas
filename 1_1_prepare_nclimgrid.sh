#!/bin/bash

dirdata="/Users/marco/Documents/dati/obs/nclimgrid"
dirwrk="./"
cd $dirdata

wget https://www.ncei.noaa.gov/thredds/fileServer/data-in-development/nclimgrid/nclimgrid_tmax.nc
wget https://www.ncei.noaa.gov/thredds/fileServer/data-in-development/nclimgrid/nclimgrid_prcp.nc

##global2california
COORD=-125,-117,34,42.5

cdo sellonlatbox,$COORD nclimgrid_tmax.nc ofile.nc
cdo selyear,1960/2021 ofile.nc nclimgrid_tmax_cal.nc
rm ofile.nc

cdo sellonlatbox,$COORD nclimgrid_prcp.nc ofile.nc
cdo selyear,1960/2021 ofile.nc nclimgrid_prcp_cal.nc
rm ofile.nc

cd $dirwrk