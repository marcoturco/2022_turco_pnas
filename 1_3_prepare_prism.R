#===============================================================================
# Description: Prepares regional averages of VPD from PRISM data
#===============================================================================

#===============================================================================
# 1). Preliminary -----
#===============================================================================

# Clean up
rm(list = ls())
graphics.off()
gc()

# Packages
# wants <- c("prism","reshape2","dplyr","ggplot2","ggmap","lubridate","stringr","stars","zoo","tidyverse","cubelyr","terra")
# install.packages('terra', repos='https://rspatial.r-universe.dev')
wants <- c("prism","terra","lubridate")
needs <- wants[!(wants %in% installed.packages()[, "Package"])]
if (length(needs))
  install.packages(needs)
lapply(wants, function(i)
  require(i, character.only = TRUE))
rm(needs, wants)

# Directories
dir <- list()
# dir$data = '/diskonfire/prism/'
dir$data = '/Users/marco/Documents/dati/obs/prism/'
dir$out = '~/Dropbox/estcena/scripts/fires_california/data_def/prism/'

#===============================================================================
# 2). get prism data
#===============================================================================

prism_set_dl_dir(dir$data)

# tmean	Mean temperature
# tmax	Maximum temperature
# tmin	Minimum temperature
# tdmean	Mean dew point temperature
# ppt	Total precipitation (rain and snow)
# vpdmin	Daily minimum vapor pressure deficit
# vpdmax	Daily maximum vapor pressure deficit

# Data is available from 1891 through present however all data prior to 1981 are downloaded as a single file.
# imonth=1
# get_prism_monthlys(type = "tmax", year = 1960:1980, mon = imonth, keepZip = FALSE)
# get_prism_monthlys(type = "tmin", year = 1960:1980, mon = imonth, keepZip = FALSE)
# get_prism_monthlys(type = "tdmean", year = 1960:1980, mon = imonth, keepZip = FALSE)

# for (iyear in 1984:2021) {
#   for (imonth in 1:12) {
#     get_prism_monthlys(type = "tmax", year = iyear, mon = imonth, keepZip = FALSE)
#     get_prism_monthlys(type = "tmin", year = iyear, mon = imonth, keepZip = FALSE)
#     get_prism_monthlys(type = "tdmean", year = iyear, mon = imonth, keepZip = FALSE)
#   }
# }


#===============================================================================
# 3). extract data over the domain
#===============================================================================

#--- starting date of the working month-year ---#
start_date <- dmy(paste0("1/1/1960"))
end_date <- dmy(paste0("1/12/2021"))
dates_ls <-substring(seq(start_date, end_date, "months"), 1, 7)
#--- remove dashes ---#
dates_prism_txt <- str_remove_all(dates_ls, "-")

var_type <- "tmax" # variable type
folder_name <- paste0("PRISM_", var_type, "_stable_4kmM3_", dates_prism_txt, "_bil") 
file_name <- paste0("PRISM_", var_type, "_stable_4kmM3_", dates_prism_txt, "_bil.bil") 
file_path <- paste0(dir$data, folder_name, "/", file_name) 
tmax_stack <- terra::rast(file_path) 
tmax_stack=terra::project(tmax_stack,"EPSG:4326")

var_type <- "tmin" # variable type
folder_name <- paste0("PRISM_", var_type, "_stable_4kmM3_", dates_prism_txt, "_bil") 
file_name <- paste0("PRISM_", var_type, "_stable_4kmM3_", dates_prism_txt, "_bil.bil") 
file_path <- paste0(dir$data, folder_name, "/", file_name) 
tmin_stack <- terra::rast(file_path) 
tmin_stack=terra::project(tmin_stack,"EPSG:4326")

var_type <- "tdmean" # variable type
folder_name <- paste0("PRISM_", var_type, "_stable_4kmM3_", dates_prism_txt, "_bil") 
file_name <- paste0("PRISM_", var_type, "_stable_4kmM3_", dates_prism_txt, "_bil.bil") 
file_path <- paste0(dir$data, folder_name, "/", file_name) 
tdmean_stack <- terra::rast(file_path) 
tdmean_stack=terra::project(tdmean_stack,"EPSG:4326")

# Latitude/Longitude WGS84 (EPSG: 4326)

data(wrld_simpl)
shp <-
  readOGR(paste0(dir$out, '/shapefiles/baileys/Baileys_ecoregions_cal.shp'))
shps = spTransform(shp,
                   CRS(
                     "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                   ))


#--- crop the entire PRISM to its domain portion---#
prism_tmax <- terra::crop(tmax_stack, shps)
prism_tmin <- terra::crop(tmin_stack, shps)
prism_tdmean <- terra::crop(tdmean_stack, shps)

# plot(prism_tmax)

sid_tmax=as.data.frame(prism_tmax, xy=TRUE)
sid_tmin=as.data.frame(prism_tmin, xy=TRUE)
sid_tdmean=as.data.frame(prism_tdmean, xy=TRUE)

lon=sort(unique(sid_tmax$x))
lat=sort(unique(sid_tmax$y))

points <- expand.grid(lon, lat)
pts = SpatialPoints(points, proj4string = CRS(proj4string(shps)))
ii <- !is.na(over(pts, shps))
inout = ii[, 1]
dim(inout) <- c(length(lon), length(lat))
inout[inout == 0] = NA
image.plot(lon, lat, inout)
plot(shps, add = TRUE)

aux_tmax = array(NA, dim = c(length(lon), length(lat), dim(tmax_stack)[3]))
aux_tmin = array(NA, dim = c(length(lon), length(lat), dim(tmax_stack)[3]))
aux_tdmean = array(NA, dim = c(length(lon), length(lat), dim(tmax_stack)[3]))


for (i in 1:length(lon)) {
  print (paste0('lon ',i,' of ',length(lon)))
  for (j in 1:length(lat)) {
    if (is.na(inout[i, j]))
      next
    for (k in 1:dim(tmax_stack)[3]) {
      iok = which(sid_tmax$x == lon[i] & sid_tmax$y == lat[j])
      if (is_empty(iok))
        next
      aux_tmax[i, j, k] = sid_tmax[iok, k + 2] 
      aux_tmin[i, j, k] = sid_tmin[iok, k + 2] 
      aux_tdmean[i, j, k] = sid_tdmean[iok, k + 2] 
    }
  }
}

image.plot(lon, lat, apply(aux_tmax,c(1,2),mean,me.rm=T))
plot(shps, add = TRUE)


# https://en.wikipedia.org/wiki/Tetens_equation
vp_tmax = 6.1078*exp(17.27*aux_tmax/(aux_tmax + 237.3)); 
vp_tmin = 6.1078*exp(17.27*aux_tmin/(aux_tmin + 237.3)); 
satvp = 6.1078*exp(17.27*aux_tdmean/(aux_tdmean + 237.3)); 
vpd=(((vp_tmax+vp_tmin)/2)-satvp)/10;

image.plot(lon, lat, apply(vpd,c(1,2),mean,me.rm=T))
plot(shps, add = TRUE)

# vp_tmax = 6.1078*exp(17.27*tmax./(tmax + 237.3)); 
# vp_tmin = 6.1078*exp(17.27*tmin./(tmin + 237.3)); 
# satvp = 6.1078*exp(17.27*tdew./(tdew + 237.3)); 
# vpd=((vp_tmax+vp_tmin)/2-satvp)/10;

tmax_prism = vector()
vpd_prism = vector()
for (im in 1:dim(vpd)[3]) {
  vpd_prism[im] = mean(vpd[, , im], na.rm = TRUE) 
  tmax_prism[im] = mean(aux_tmax[, , im], na.rm = TRUE) 
}

plot.ts(as.vector(vpd_prism))
save(vpd_prism, file = file.path(dir$out,"vpd_prism.RData"))
writeMat(file.path(dir$out, paste("vpd_prism.mat",sep="")), vpd_prism=vpd_prism)

plot.ts(as.vector(tmax_prism))
save(tmax_prism, file = file.path(dir$out,"tmax_prism.RData"))
writeMat(file.path(dir$out, paste("tmax_prism.mat",sep="")), tmax_prism=tmax_prism)

