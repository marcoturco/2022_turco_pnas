#===============================================================================
# Description: Prepares region-level nclimgrid data, regional averages of tmax 
# and prec
#===============================================================================

#===============================================================================
# 1). Preliminary -----
#===============================================================================

# Clean up
rm(list = ls())
graphics.off()
gc()

# Packages
wants <- c("raster","rgdal", "maptools", "fields", "maps", "ncdf4",
           "rgeos","R.matlab")
needs <- wants[!(wants %in% installed.packages()[, "Package"])]
if (length(needs))
  install.packages(needs)
lapply(wants, function(i)
  require(i, character.only = TRUE))
rm(needs, wants)

# Directories
dir <- list()
dir$root = '~/Dropbox/estcena/scripts/fires_california/'
dir$data = paste0(dir$root,'/data_def/')
dir$out = paste0(dir$data,'/nclimgrid/')
dir$obs = '~/Documents/dati/obs/nclimgrid/'

# Misc
misc <- list()

#===============================================================================
# 2). shapefile -----
#===============================================================================

data(wrld_simpl)
shp <-
  readOGR(paste0(dir$data, '/shapefiles/baileys/Baileys_ecoregions_cal.shp'))
shp = spTransform(shp,
                  CRS(
                    "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                  ))

#===============================================================================
# 3). load the nclimgrid data -----
#===============================================================================
data(wrld_simpl)
#tmax
fname <- file.path(dir$obs, 'nclimgrid_tmax_cal.nc')
obs.nc <- nc_open(fname)
tmax <- ncvar_get(obs.nc, "tmax")
tmax[tmax == "NaN"] <- NA
obs.nc$dim$lon$vals -> lon
obs.nc$dim$lat$vals -> lat
lat = rev(lat)
tmax = tmax[, ncol(tmax):1, ]
image.plot(lon, lat, apply(tmax, c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)
points <- expand.grid(lon, lat)
pts = SpatialPoints(points, proj4string = CRS(proj4string(shp)))



#pr
fname <- file.path(dir$obs, 'nclimgrid_prcp_cal.nc')
obs.nc <- nc_open(fname)
prcp <- ncvar_get(obs.nc, "prcp")
prcp[prcp == "NaN"] <- NA
prcp = prcp[, ncol(prcp):1, ]
image.plot(lon, lat, apply(prcp, c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

#spatial averages
tx = matrix(data = NA,
            nrow = dim(tmax)[3],
            ncol = dim(shp)[1])
pr = matrix(data = NA,
            nrow = dim(prcp)[3],
            ncol = dim(shp)[1])

for (ireg in 1:dim(shp)[1]) {
  aux1 = shp[ireg,]@polygons[[1]]@Polygons[[1]]@coords
  inout = point.in.polygon(points[, 1], points[, 2], aux1[, 1], aux1[, 2])
  dim(inout) <- c(nrow(tmax), ncol(tmax))
  inout[inout == 0] = NA
  if (sum(inout, na.rm = TRUE) > 0) {
    image.plot(lon, lat, inout)
    tmax2 = tmax
    prcp2 = prcp
    for (i in 1:dim(tmax)[3]) {
      tmax2[, , i] = tmax[, , i] * inout
      prcp2[, , i] = prcp[, , i] * inout
      tx[i, ireg] = mean(tmax2[, , i], na.rm = TRUE)
      pr[i, ireg] = mean(prcp2[, , i], na.rm = TRUE)
    }
  }
}

nclimgrid = array(data = NA, dim = c(2, dim(tx)[1])) 
nclimgrid[1, ] = tx
nclimgrid[2, ] = pr

save(nclimgrid, file = file.path(dir$out, paste0("nclimgrid.RData")))
writeMat(file.path(dir$out, paste0("nclimgrid.mat")), nclimgrid = nclimgrid)