#### Create Basins ####
## This file read in the matlab file and create polygon shapefiles for the basins


## First load the matlab file
library(R.matlab)
basinmat <- readMat("Z:/WP3-Ocean/Data/BasinMasks/regionMap.mat")
image(t(basinmat$reg), asp = 0.8)
keys <- c("SAt", "ISPac", "EPac", "STNAt", "SPNAt", "NWPac")
fullkeys <- sapply(sapply(basinmat$key, "["), "[")

## Create a raster from the matrix
library(raster)
basin_raster <- flip(raster(basinmat$reg, xmn=0, xmx=360, ymn= -90, ymx = 90), 2)

basinsPolys <- lapply(1:6, function(x) {rasterToPolygons(basin_raster, fun = function(k) {x==k}, dissolve = TRUE)} )

Plist <- list()
for (i in 1:6){
plist <- basinsPolys[[i]]@polygons[[1]]@Polygons
Plist[[i]] <- Polygons(plist, ID = keys[i])
}

allbasins <- SpatialPolygons(Plist)
plot(allbasins, asp = 1, col = 1:6)
allbasinsdf <- SpatialPolygonsDataFrame(allbasins, data = data.frame(names = fullkeys), match.ID = FALSE)
## Save this as a shapefile
library(rgdal)
writeOGR(allbasinsdf, dsn = "Z:/WP3-Ocean/Data/BasinMasks/basinPolys", layer = "Basins", driver = "ESRI Shapefile")
save(allbasinsdf, file ="Z:/WP1-BHM/maps/Ocean/basins.rda")
