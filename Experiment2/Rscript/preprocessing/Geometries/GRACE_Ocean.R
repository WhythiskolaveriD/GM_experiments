#### Create GRACE Ocean Polygons ####
## This file create the Ocean polygons based on MS's superMascon data
library(rgdal)
library(rgeos)
## Load MS superMascon definition
gracetab <- read.table(file = "Z:/WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_loc_v03.txt", header = FALSE, skip = 1)
## - V1: id
## - V2: area
## - V3: long_center
## - V4: lat_center
## - V5 - V78: long/lat

## Create Polygons info dataframe
datainfo <- gracetab[,1:4]
names(datainfo) <- c("ID", "area", "longC", "latC")

## Create Polygons based on longlat points
datalon <- gracetab[,c(rep(FALSE, 4), rep(c(TRUE,FALSE), 37))]
datalat <- gracetab[,c(rep(FALSE, 4), rep(c(FALSE,TRUE), 37))]
LonLat <- lapply(1:4278, function(x) na.omit(cbind(as.numeric(datalon[x,]), as.numeric(datalat[x,]))))

srl <- lapply(1:4278, function(x) Polygon(coords =  LonLat[[x]]))
Srl <- lapply(1:4278, function(x)Polygons(srl[x], ID = datainfo$ID[x]))
SP_Supermas <- SpatialPolygons(Srl, proj4string = CRS("+proj=longlat"))
plot(SP_Supermas)

## Next separate Ocean and Land Masscons
## Load the land Polygons
land0 <- readOGR(dsn = "Z:/WP1-BHM/maps/ne_110m_land", layer = "ne_110m_land")
## Remove small islands in the low lattidure(lat < 70) regrions (area < 300*300km)
areas <- sapply(land0@polygons, slot, "area")
lats <- sapply(land0@polygons, function(x){slot(x, "labpt")[2]})
smallareas <- areas <= 10
lowlats <- lats <= 70
lowsmall <- lowlats & smallareas
newids <- which(!lowsmall)
Land <- SpatialPolygons(land0@polygons[newids], proj4string = CRS(proj4string(land0)))

## Find Supermascons on land
mas_lonlat <- cbind(datainfo$longC, datainfo$latC)
mas_lonlat[,1] <- ifelse(mas_lonlat[,1] > 180, mas_lonlat[,1]-360, mas_lonlat[,1])
masC_sp <- SpatialPoints(coords = mas_lonlat, proj4string = CRS(proj4string(Land)))
overid <- over(masC_sp, Land)
oceanid <- which(is.na(overid))
landid <-  which(!is.na(overid))
plot(masC_sp[oceanid, ], col = "blue")
plot(masC_sp[landid, ], col = "orange", add = T)

## select the ocean polygons
Ocean_poly <- SP_Supermas[oceanid,]
Land_poly <- SP_Supermas[landid,]
plot(Ocean_poly, col = "blue")
plot(Land_poly, col = "orange", add = T)

Ocean_spf <- SpatialPolygonsDataFrame(Ocean_poly, datainfo[oceanid,])
Land_spf <- SpatialPolygonsDataFrame(Land_poly, datainfo[landid,])

## Save them as shape files and rda files
writeOGR(Ocean_spf, dsn = "Z:/WP2-SolidEarth/Bramha/Datsets/catchment_shp", layer = "OceanSupMas", driver = "ESRI Shapefile")
save(Ocean_spf, file ="Z:/WP1-BHM/maps/GRACE/oceanMas_poly.rda")

writeOGR(Land_spf, dsn = "Z:/WP2-SolidEarth/Bramha/Datsets/catchment_shp", layer = "LandSupMas", driver = "ESRI Shapefile")
save(Land_spf, file ="Z:/WP1-BHM/maps/GRACE/landMas_poly.rda")


## We also want to remove the masscons near the land boundaries since the signals are contaminated
## Therefoer we create land bufferzones

### Create buffer zones for large land polygons.
## -- The land polygons are devided into large and small groups.
## -- Assuming the correlation lenght on the land can be as large as 600km, 
## -- then "Small" is defined by area < 400,000 km^2. This also cover most of small islands.
## -- The rest are "large", and we remove the polygon hole that represent the Caspian Sea.
### Buffer zones is 200km out of land

### Remove the Caspian Sea
holes <- which(sapply(Land@polygons, function(x) length(x@Polygons) > 1))
caspian <- no_caspian <- Land@polygons[[holes]]
no_caspian@Polygons <- no_caspian@Polygons[1]
no_caspian@plotOrder <- 1L

## The land SpatialPolygons with no Caspian sea
Land2 <- Land
Land2@polygons[[holes]] <- no_caspian 
slot(Land2, "polygons") <- lapply(slot(Land2, "polygons"), "comment<-", NULL) 

## id of Antarcitca
plot(Land2[1:3,])


### Create the buffer for large polygons
## All the polygons except Antarctica
p1 <- Land2[-c(1:3),]
## Transform the CRS of the large polygons for the bufferring operation
## Use a proper CRS: epsg:3857, see https://epsg.io/3857
## Works for lats within 85 -- need another one for Antarctica
p1 <- spTransform(p1, CRS("+init=epsg:3857"))
bp1 <- gBuffer(p1, width = 200*1000 )
## Remove small area regions and holes in the buffers
holesid <- which(sapply(bp1@polygons[[1]]@Polygons, slot, "hole"))
bp1 <- SpatialPolygons(list(Polygons(bp1@polygons[[1]]@Polygons[-holesid], ID = 1)),proj4string = CRS("+init=epsg:3857") )

## Use another one for Antarctica, see https://epsg.io/5481
p2 <- Land2[1:3,]
p2 <- spTransform(p2, CRS("+init=epsg:5481"))
bp2 <- gBuffer(p2, width = 200*1000 )
## !!!Note: no need to project them back! Project the points when using over!!!

## Find Supermascons on land
mas_lonlat <- cbind(datainfo$longC, datainfo$latC)
mas_lonlat[,1] <- ifelse(mas_lonlat[,1] > 180, mas_lonlat[,1]-360, mas_lonlat[,1])
masC_sp <- SpatialPoints(coords = mas_lonlat, proj4string = CRS(proj4string(Land)))
masC_sp1 <- spTransform(masC_sp, CRS(proj4string(bp1)))
masC_sp2 <- spTransform(masC_sp, CRS(proj4string(bp2)))
overid1 <- over(masC_sp1, bp1)
overid2 <- over(masC_sp2, bp2)
idall <- 1:length(masC_sp)
landid <-  c(which(!is.na(overid1)), which(!is.na(overid2)))
oceanid <- idall[-landid]

plot(masC_sp[oceanid, ], col = "blue")
plot(masC_sp[landid, ], col = "orange", add = T)

## select the ocean polygons
Ocean_poly <- SP_Supermas[oceanid,]
Land_poly <- SP_Supermas[landid,]
plot(Ocean_poly, col = "blue")
plot(Land_poly, col = "orange", add = T)

Ocean_spf <- SpatialPolygonsDataFrame(Ocean_poly, datainfo[oceanid,])
Land_spf <- SpatialPolygonsDataFrame(Land_poly, datainfo[landid,])

## Save them as shape files and rda files
writeOGR(Ocean_spf, dsn = "Z:/WP2-SolidEarth/Bramha/Datsets/catchment_shp", layer = "OceanSupMas", driver = "ESRI Shapefile")
save(Ocean_spf, file ="Z:/WP1-BHM/maps/GRACE/oceanMas_poly.rda")

writeOGR(Land_spf, dsn = "Z:/WP2-SolidEarth/Bramha/Datsets/catchment_shp", layer = "LandSupMas", driver = "ESRI Shapefile")
save(Land_spf, file ="Z:/WP1-BHM/maps/GRACE/landMas_poly.rda")

