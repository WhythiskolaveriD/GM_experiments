#### Create GRACE Catchment Regions ####
## This file read in the matlab file and create polygon shapefiles for the grace cathcment regions

## First load the matlab file
library(R.matlab)
grace_mat <- readMat("Z:/WP2-SolidEarth/Bramha/Datsets/catchment_polygons.mat")
id_useful <- unlist(readMat("Z:/WP2-SolidEarth/Bramha/Datsets/index_useful_catchment.mat"))
## Looking at the original mat file in matlab, the file contains the following 9 fields of 405 entries
## The fields names are 1.Geometry, 2.BoundingBox, 3.Lon, 4.Lat, 5.Catch_id, 6.name, 7.rasterarea, 8.DB_ID, 9.Lattitude
## Create the the index for useful fields
BBox_idx <- which((1:3645 %% 9) == 2)
Lon_idx <- which((1:3645 %% 9) == 3)
Lat_idx <- which((1:3645 %% 9) == 4)
ID_idx <- which((1:3645 %% 9) == 5)
name_idx <- which((1:3645 %% 9) == 6)
area_idx <- which((1:3645 %% 9) == 7)
DBID_idx <- which((1:3645 %% 9) == 8)
LAT_idx <- which((1:3645 %% 9) == 0)

## Next We use the relavent information to create spatial polygons for each catchment regions and then the dataframe
## Create the polygons
library(sp)
IDs <- unlist(grace_mat$basin[ID_idx])
Lons <- grace_mat$basin[Lon_idx]
Lats <- grace_mat$basin[Lat_idx]
NA_num <- sapply(Lons, function(x) sum(is.na(x)))
natab <- table(NA_num)

## 1 NA points -- remove the NAs
Lonlat1 <- lapply(which(NA_num == 1), function(x) na.omit(cbind(t(Lons[[x]]), t(Lats[[x]]))))

## More than 1 NA points -- check the sequence length
seql <- list()
ids2 <- which(NA_num > 1)
for (i in 1:length(ids2)){
  naid <- which(is.na(Lons[[ids2[i]]]))
  naid0 <- c(0, naid[1:(length(naid) - 1)])
  seqli <- naid - naid0
  seql[[i]] <- cbind(naid0, naid, seqli)
}

## Keep the longest sequence
LonLat2 <- lapply(1:length(ids2), 
                  function(x) {
                    lons <- Lons[[ids2[x]]]
                    lats <- Lats[[ids2[x]]]
                    seqtemp <- seql[[x]]
                    maxid <- which(seqtemp[,3] == max(seqtemp[,3]))
                    id0 <- seqtemp[maxid,1]+1
                    id1 <- seqtemp[maxid,2]-1
                    coord <- cbind(lons[id0:id1], lats[id0:id1])
                    return(coord)})
## Assemble the two
Lonlat <- list()
Lonlat[which(NA_num == 1)] <- Lonlat1
Lonlat[ids2] <- LonLat2

## Build the SpatialPolygon
srl <- lapply(1:405, function(x) Polygon(coords =  Lonlat[[x]]))
Srl <- lapply(1:405, function(x)Polygons(srl[x], ID = IDs[x]))
SP1 <- SpatialPolygons(Srl, proj4string = CRS("+proj=longlat"))
plot(SP1)

## Build the Dataframe associate with the polygons
datainfo <- data.frame(name = unlist(grace_mat$basin[name_idx]),
                       area = unlist(grace_mat$basin[area_idx]),
                       DBID = unlist(grace_mat$basin[DBID_idx]),
                       LAT  = unlist(grace_mat$basin[LAT_idx]))
row.names(datainfo) <- as.character(unlist(grace_mat$basin[ID_idx]))
SPdf <- SpatialPolygonsDataFrame(SP1, datainfo)
## Save these as shapefile for future use
library(rgdal)
writeOGR(SPdf, dsn = "Z:/WP2-SolidEarth/Bramha/Datsets/catchment_shp", layer = "landcatchment", driver = "ESRI Shapefile")

## Also save the SpatialPolygonDataFrame in an rda file
save(SPdf, file ="Z:/WP1-BHM/maps/GRACE/landcatchment_all.rda")

## Next we select the useful regions
SPdf_useful <- SPdf[id_useful,]
plot(SPdf_useful)
## Save these as shapefile and rda file
writeOGR(SPdf_useful, dsn = "Z:/WP2-SolidEarth/Bramha/Datsets/catchment_shp", layer = "landcatchment_usefull", driver = "ESRI Shapefile")
save(SPdf_useful, file ="Z:/WP1-BHM/maps/GRACE/landcatchment_usefull.rda")

## Define Basins for Antarctica
Antarctica <- readOGR(dsn = "Z:/WP1-BHM/maps/Antarctica/ANT_Basins_IMBIE2_v1.6", layer = "ANT_Basins_IMBIE2_v1.6")
## Remove islands
Antarctica2 <- Antarctica[2:19,]

## Define basins for Greenland
Greenland <- readOGR(dsn = "Z:/WP1-BHM/maps/Greenland/GRE_Basins_IMBIE2_v1.3", layer = "GRE_Basins_IMBIE2_v1.3")
## Remove small islands
Greenland2 <- Greenland[52:57,]
Greenland2@polygons[[5]]@Polygons <- Greenland2@polygons[[5]]@Polygons[1]

save(SPdf_useful, Antarctica2, Greenland2, file ="Z:/WP1-BHM/maps/GRACE/landcatchment_usefull.rda")
