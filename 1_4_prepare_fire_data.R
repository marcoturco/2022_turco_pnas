#===============================================================================
# Description: Prepares region-level fire data
#===============================================================================

#===============================================================================
# 1). Preliminary -----
#===============================================================================

# Clean up
rm(list = ls())
graphics.off()
gc()

# Packages
wants <- c("raster","rgdal", "rgeos","R.matlab")
needs <- wants[!(wants %in% installed.packages()[, "Package"])]
if (length(needs))
  install.packages(needs)
lapply(wants, function(i)
  require(i, character.only = TRUE))
rm(needs, wants)

# Directories
dir <- list()
dir$data = '~/Dropbox/estcena/scripts/fires_california/data_def/'
dir$fire = '~/Documents/dati/fire_us/'
dir$out = paste0(dir$data,'/nclimgrid/')

# Misc
misc <- list()
misc$years_mtbs <- 1984:2021
misc$years_frap <- 1971:2021

#===============================================================================
# 2). MTBS data -----
#===============================================================================

# The input file geodatabase
mtbs <- readOGR(paste0(dir$fire,"/mtbs21_sierra_northcoast.shp"), verbose = FALSE)
# summary(mtbs)
mtbs = spTransform(mtbs,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

BA = array(data = 0, dim = c(1, length(misc$years_mtbs) * 12))

anni=as.numeric(format(as.Date(mtbs@data$Ig_Date),"%Y"))
mesi=as.numeric(format(as.Date(mtbs@data$Ig_Date),"%m"))


k = 0
for (iyear in misc$years_mtbs) {
  for (imonth in 1:12) {
    k = k + 1
    idx2 = which(anni == iyear &
                   mesi == imonth)
    if (length(idx2) >= 1) {
      
      
      for (iidx2 in 1:length(idx2)) {
        a=unlist(mtbs[idx2[iidx2],]@polygons)
        BA[1, k] = sum(BA[1, k],a[[1]]@Polygons[[1]]@area,na.rm = TRUE)
      } 
    }
  }
}

plot.ts(as.vector(BA))
save(BA, file = file.path(dir$data,"fires/mtbs21_sierra_northcoast.RData"))
writeMat(file.path(dir$data, paste("fires/mtbs21_sierra_northcoast.mat",sep="")), BA=BA)

# Over northen california
# The input file  (previously cut in QGIS)
frap <- readOGR(paste0(dir$fire,"fire21_sierra_northcoast.shp"), verbose = FALSE)
frap = spTransform(frap,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

BA = array(data = 0, dim = c(1, length(misc$years_frap) * 12))
anni=as.numeric(format(as.Date(frap@data$ALARM_D),"%Y"))
mesi=as.numeric(format(as.Date(frap@data$ALARM_D),"%m"))
giorni=as.numeric(format(as.Date(frap@data$ALARM_D),"%d"))

fires=frap@data$GIS_ACR

k = 0
for (iyear in misc$years_frap) {
  for (imonth in 1:12) {
    k = k + 1
    idx2 = which(anni == iyear &
                   mesi == imonth)
    if (length(idx2) >= 1) {
      BA[1, k] = sum(frap[idx2,]@data$GIS_ACR)
    }
  }
}

writeMat(file.path(dir$data, paste("fires/all_fire21_sierra_northcoast.mat",sep="")), BA=BA)
