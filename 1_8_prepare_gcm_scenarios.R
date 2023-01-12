#===============================================================================
# Description: Prepares regional averages of GCM tasmax
#===============================================================================

#===============================================================================
# 1). Preliminary -----
#===============================================================================

# Clean up
rm(list = ls())
graphics.off()
gc()

# Packages
wants <- c("rgdal",
           "maptools",
           "fields",
           "maps",
           "ncdf4",
           "rgeos",
           "R.matlab")
needs <- wants[!(wants %in% installed.packages()[, "Package"])]
if (length(needs))
  install.packages(needs)
lapply(wants, function(i)
  require(i, character.only = TRUE))
rm(needs, wants)

# Directories
dir <- list()
# dir$root <- dirname(getwd())
dir$root = '~/Dropbox/estcena/scripts/fires_california/'
dir$data = '~/Dropbox/estcena/scripts/fires_california/data_def/'
dir$gcm = '~/Documents/dati/gcm_data_cal/'
dir$out = '~/Dropbox/estcena/scripts/fires_california/data_def/gcms/' #tasmax/ssp245/'

# Misc
misc <- list()
misc$years_hist <- 1950:2014 # years hist-nat simulations
misc$years_fut <- 2015:2100 # years future simulations

#===============================================================================
# 2). shapefile -----
#===============================================================================

data(wrld_simpl)
shp <-
  readOGR(paste0(dir$data, '/shapefiles/baileys/Baileys_ecoregions_cal.shp'))
shps = spTransform(shp,
                   CRS(
                     "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                   ))

#===============================================================================
# 3). load the projections data -----
#===============================================================================

listGCMs = list.files(path = paste0(dir$gcm, 'gcm_1_x_1/tasmax/'),
                      pattern = "historical(.*)nc$")

length(listGCMs) #453

listGCMs[1]
# extract lat lon and mask
filename=paste0(dir$gcm, 'gcm_1_x_1/tasmax/', listGCMs[1])
gcm.nc <- nc_open(filename)
gcm.nc$dim$lon$vals -> lon
gcm.nc$dim$lat$vals -> lat
ni = length(lon)
nj = length(lat)
points <- expand.grid(lon, lat)
pts = SpatialPoints(points, proj4string = CRS(proj4string(shps)))
ii <- !is.na(over(pts, shps))
inout = ii[, 1]
dim(inout) <- c(ni, nj)
inout[inout == 0] = NA
image.plot(lon, lat, inout)
plot(shps, add = TRUE)

# tasmax
variab = "tasmax"
kmod = 0
model_name=character(length(listGCMs)) 
historical <- matrix(data=NA,nrow=length(misc$years_hist)*12,ncol=length(listGCMs))
ssp245 <- matrix(data=NA,nrow=length(misc$years_fut)*12,ncol=length(listGCMs))
ssp585=ssp245

for (ifile in 1:length(listGCMs)) {
  ## check if historical, ss245 and ssp585 simulation exist for the same member of the same model
  # "tasmax_CMIP6Amon_ACCESS-CM2_historical_r1i1p1f1_mon_1950-2014_CA_Yizhou.nc"
  # tasmax_CMIP6Amon_UKESM1-0-LL_ssp585_r8i1p1f2_mon_2015-2100_CA_Yizhou.nc
  namefile_hist = paste0(dir$gcm, 'gcm_1_x_1/tasmax/', listGCMs[ifile])
  aux = gsub('1950-2014', '2015-2100', listGCMs[ifile])
  namefile_ssp245 = paste0(dir$gcm,
                           'gcm_1_x_1/tasmax/',
                           gsub('historical', 'ssp245', aux))
  namefile_ssp585 = paste0(dir$gcm,
                           'gcm_1_x_1/tasmax/',
                           gsub('historical', 'ssp585', aux))
  
  
  if (file.exists(namefile_hist) &&
      file.exists(namefile_ssp245) &&
      file.exists(namefile_ssp585)) {
    kmod = kmod + 1
    aux = gsub('tasmax_CMIP6Amon_', '', aux)
    aux = gsub('_historical_', '_', aux)
    aux = gsub('_mon_2015-2100_CA_Yizhou.nc', '', aux)
    model_name[kmod]=aux
    #historical
    gcm.nc <- nc_open(namefile_hist)
    gcm <- ncvar_get(gcm.nc, variab)
    gcm[gcm == "1.00000002004088e+20"] <- NA
    image.plot(lon, lat, apply(gcm, c(1, 2), mean, na.rm = TRUE))
    plot(wrld_simpl, add = TRUE)
    plot(shps, add = TRUE)
    gcm2 = gcm
    tx_gcm = vector()
    for (im in 1:dim(gcm)[3]) {
      gcm2[, , im] = inout * gcm[, , im]
      tx_gcm[im] = mean(gcm2[, , im], na.rm = TRUE)
    }
    historical[,kmod]=tx_gcm
    #ssp245
    gcm.nc <- nc_open(namefile_ssp245)
    gcm <- ncvar_get(gcm.nc, variab)
    gcm[gcm == "1.00000002004088e+20"] <- NA
    image.plot(lon, lat, apply(gcm, c(1, 2), mean, na.rm = TRUE))
    plot(wrld_simpl, add = TRUE)
    plot(shps, add = TRUE)
    gcm2 = gcm
    tx_gcm = vector()
    for (im in 1:dim(gcm)[3]) {
      gcm2[, , im] = inout * gcm[, , im]
      tx_gcm[im] = mean(gcm2[, , im], na.rm = TRUE)
    }
    ssp245[,kmod]=tx_gcm
    #ssp585
    gcm.nc <- nc_open(namefile_ssp585)
    gcm <- ncvar_get(gcm.nc, variab)
    gcm[gcm == "1.00000002004088e+20"] <- NA
    image.plot(lon, lat, apply(gcm, c(1, 2), mean, na.rm = TRUE))
    plot(wrld_simpl, add = TRUE)
    plot(shps, add = TRUE)
    gcm2 = gcm
    tx_gcm = vector()
    for (im in 1:dim(gcm)[3]) {
      gcm2[, , im] = inout * gcm[, , im]
      tx_gcm[im] = mean(gcm2[, , im], na.rm = TRUE)
    }
    ssp585[,kmod]=tx_gcm
  }
}

namefile = paste0(dir$data,'gcms/list_gcms_future.mat');
writeMat(namefile, model_name  = model_name)

namefile = paste0(dir$data,'gcms/tasmax-historical-1950-2014-monthly.mat');
writeMat(namefile, historical = historical)

namefile = paste0(dir$data,'gcms/tasmax-ssp245-2015-2100-monthly.mat');
writeMat(namefile, ssp245 = ssp245)

namefile = paste0(dir$data,'gcms/tasmax-ssp585-2015-2100-monthly.mat');
writeMat(namefile, ssp585 = ssp585)