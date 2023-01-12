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
                      pattern = "(.*)nc$")

# tasmax
variab = "tasmax"
for (ifile in 1:length(listGCMs)) {
  ## load ctrl
  # tasmax_CMIP6Amon_UKESM1-0-LL_historical_r9i1p1f2_mon_1950-2014_CA.nc
  namefile=paste0(dir$gcm, 'gcm_1_x_1/tasmax/',listGCMs[ifile])
  gcm.nc <- nc_open(namefile)
  gcm <- ncvar_get(gcm.nc, variab)
  gcm[gcm == "1.00000002004088e+20"] <- NA
  gcm.nc$dim$lon$vals -> lon
  gcm.nc$dim$lat$vals -> lat
  image.plot(lon, lat, apply(gcm, c(1, 2), mean, na.rm = TRUE))
  plot(wrld_simpl, add = TRUE)
  plot(shps, add = TRUE)
  
  ni = dim(gcm)[1]
  nj = dim(gcm)[2]
  points <- expand.grid(lon, lat)
  pts = SpatialPoints(points, proj4string = CRS(proj4string(shps)))
  ii <- !is.na(over(pts, shps))
  inout = ii[, 1]
  dim(inout) <- c(nrow(gcm), ncol(gcm))
  inout[inout == 0] = NA
  image.plot(lon, lat, inout)
  plot(shps, add = TRUE)
  
  gcm2 = gcm
  tx_gcm = vector()
  for (im in 1:dim(gcm)[3]) {
    gcm2[, , im] = inout * gcm[, , im]
    tx_gcm[im] = mean(gcm2[, , im], na.rm = TRUE) 
  }
  
  namefile2 = paste0(substr(namefile,1,nchar(namefile)-23),'.RData') 
  # print(namefile)
  # save(tx_gcm, file = namefile2)
  writeMat(paste0(gsub(".RData", "", namefile2), '.mat'), tx_gcm = tx_gcm)
}


## pr
listGCMs = list.files(path = paste0(dir$gcm, 'gcm_1_x_1/pr/'),
                      pattern = "(.*)nc$")

# pr
variab = "pr"
for (ifile in 1:length(listGCMs)) {
  print(paste0('file ',ifile,' of ',length(listGCMs)))
  ## load ctrl
  # tasmax_CMIP6Amon_UKESM1-0-LL_historical_r9i1p1f2_mon_1950-2014_CA.nc
  namefile=paste0(dir$gcm, 'gcm_1_x_1/pr/',listGCMs[ifile])
  gcm.nc <- nc_open(namefile)
  gcm <- ncvar_get(gcm.nc, variab)
  # gcm[gcm == "1.00000002004088e+20"] <- NA
  gcm.nc$dim$lon$vals -> lon
  gcm.nc$dim$lat$vals -> lat
  image.plot(lon, lat, apply(gcm, c(1, 2), mean, na.rm = TRUE))
  plot(wrld_simpl, add = TRUE)
  plot(shps, add = TRUE)
  
  ni = dim(gcm)[1]
  nj = dim(gcm)[2]
  points <- expand.grid(lon, lat)
  pts = SpatialPoints(points, proj4string = CRS(proj4string(shps)))
  ii <- !is.na(over(pts, shps))
  inout = ii[, 1]
  dim(inout) <- c(nrow(gcm), ncol(gcm))
  inout[inout == 0] = NA
  image.plot(lon, lat, inout)
  plot(shps, add = TRUE)
  
  gcm2 = gcm
  pr_gcm = vector()
  for (im in 1:dim(gcm)[3]) {
    gcm2[, , im] = inout * gcm[, , im]
    pr_gcm[im] = mean(gcm2[, , im], na.rm = TRUE) 
  }
  
  namefile2 = paste0(substr(namefile,1,nchar(namefile)-23),'.RData') 
  # namefile2 = paste0(substr(namefile,1,nchar(namefile)-2),'RData')
  # print(namefile)
  # save(pr_gcm, file = namefile2)
  writeMat(paste0(gsub(".RData", "", namefile2), '.mat'), pr_gcm = pr_gcm)
}
