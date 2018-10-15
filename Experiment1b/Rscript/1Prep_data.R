###################################################################
##                  Data Processing and Mesh                     ##
## This script read in and process the data and produce the mesh ##
###################################################################
#### 0 load libraries
library(GEOmap)
library(INLA)
if(Sys.info()['sysname'] == "Linux"){INLA:::inla.dynload.workaround()}

#### 1 Read in GIA data
###################################################################
GIA_prior <- read.table(GIA_file, header = T)

#### 2 Genreate the mesh and assign the mean value
###################################################################
## Suppose we want predict at a 1 degree resolution ~ 100km
## In this experiment the residuals should be almost independent at 100km apart...
Mesh_GIA <- inla.mesh.create(globe= 64)

## assign the polygon value to the mesh vertices
polycoords <- GIA_prior[,c(6:13, 6,7)]
plist <- lapply(GIA_prior$ID, function(x) Polygons(list(Polygon(cbind(lon = as.numeric(polycoords[x, c(1,3,5,7,9)]), 
                                                                      lat = as.numeric(polycoords[x, c(2,4,6,8,10)])))), ID = x))
Plist <- SpatialPolygons(plist, proj4string = CRS("+proj=longlat"))
MlocLL <- Lxyz2ll(list(x=Mesh_GIA$loc[,1], y = Mesh_GIA$loc[,2], z = Mesh_GIA$loc[,3]))
MlocLL$lon <- ifelse(MlocLL$lon <0, MlocLL$lon + 360, MlocLL$lon)
MlocLL$lon <- ifelse(MlocLL$lon > 359.5, MlocLL$lon - 360, MlocLL$lon)
M_sp <- SpatialPoints(data.frame(lon = MlocLL$lon, lat = MlocLL$lat), proj4string = CRS("+proj=longlat")) #This convert GIA_ice6g a SpatialPointDataFrame
Midx <- over(M_sp, Plist)

GIA_mu <- GIA_prior$trend[Midx]

#### 3 Read GPS data
###################################################################
GPS_obs <- read.table(GPS_file, header = T)
## Transform the GPS location to 3D unit ball
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$lat, lon = GPS_obs$lon))

## Find the prior mean at the GPS locations
GPSY <- GPS_obs$lat
GPSx <- ifelse(GPS_obs$lon < 0, GPS_obs$lon + 360, GPS_obs$lon)
GPSxx <- ifelse(GPSx > 359.5, GPSx - 360, GPSx)
GPS_sp <- SpatialPoints(data.frame(lon = GPSxx, lat = GPSY), proj4string = CRS("+proj=longlat"))
Midx2 <- over(GPS_sp, Plist)
GPS_mu <- GIA_prior$trend[Midx2]

## Detrend the mean to get the GPS_data
GPS_data <- GPS_obs$trend - GPS_mu

#### 4 Save all initial built up objects
###################################################################
#save(Mesh_GIA, GIA_mu, GPS_loc, GPS_data, Plist, file = paste0(outdir, GIA_input, "mesh.RData"))

