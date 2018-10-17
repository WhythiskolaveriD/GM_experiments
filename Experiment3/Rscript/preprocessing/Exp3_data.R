###### This script load and process the data required by Experiment 3
####################################################################################################
datamap <- paste0(dd_save, expname, "_datamap.rda")

if(file.exists(datamap)){
  load(datamap)
}else{
  #### Load raw data
  ## The North America GPS data
  GPS <- read.table(paste0(dd_root,"WP2-SolidEarth/BHMinputs/GPS/GPS_v04b_NorthAmerica.txt"), header = T)
  coordinates(GPS) <- c("lon", "lat")
  ## The GSFC GRACE masscon data
  grace_nc <- nc_open(paste0(dd_root, "WP2-SolidEarth/Bramha/BHMinput/GRACE_VCDtrends_JPLRL06.nc"))
  print(grace_nc)
  lon <- ncvar_get(grace_nc, "Long")
  lat <- ncvar_get(grace_nc, "Lat")
  trend_grace <- ncvar_get(grace_nc, "trend") #note that there are NAs for land datat
  err_grace <- ncvar_get(grace_nc, "uncert")
  
  grace_raw <- data.frame(trend = as.numeric(trend_grace), std = as.numeric(err_grace), lon= lon, lat = lat)
  coordinates(grace_raw) <- c("lon", "lat")
  
  ## The ice6g GIA solution
  gia_ewh <- read.table(paste0(dd_root, "WP2-SolidEarth/BHMinputs/GIA/gia_ICE6G_vcd.txt"), sep = "", skip = 1)
  names(gia_ewh) <- c("lon", "lat", "trend")
  gia_ewh$lon <- ifelse(gia_ewh$lon < 0, gia_ewh$lon + 360, gia_ewh$lon)
  gia_ewh$error <- abs(0.2*gia_ewh$trend) # assume the errors are 20% of the gia trend
  coordinates(gia_ewh) <- c("lon", "lat")
  
  
  ####################################################################################################
  
  
  #### Make the North American SpatialPolygon
  Continents <- readOGR(dsn = paste0(dd_root, "WP1-BHM/maps/Continents"), layer = "continent")
  NAm <- as( Continents[2,], "SpatialPolygons")
  
  ## Remove Greenland a part of small islands.
  summary(NAm)
  areas <- sapply(NAm@polygons[[1]]@Polygons, slot, "area")
  greenid <- order(areas, decreasing = TRUE)[2]
  smallid <- which(areas < 5)
  NorthA <- NAm@polygons[[1]]@Polygons[-c(greenid,smallid)]
  polycent <- sapply(NorthA, slot, "labpt")
  
  NorthA <- lapply(1:16, function(x) Polygons(list(NorthA[[x]]), ID = x))
  NorthA <- SpatialPolygons(NorthA, proj4string = CRS(proj4string(NAm)))
  
  ## Remove north Canadian islands and small island in the south
  NC <- which(polycent[2, ] > 70)
  SA <- which(polycent[2, ] < 25)
  cols <- rep(NA, length(polycent))
  ## plot
  #cols[NC] <- 2
  #cols[SA] <- 3
  #plot(NorthA, col = cols)
  
  NorthA <- NorthA[-c(NC, SA)]
  
  ## Transform to the "North America Albers Equal Area Conic projection" 
  ## see http://spatialreference.org/ref/esri/north-america-albers-equal-area-conic/
  NorthA_t <- spTransform(NorthA, CRS("+init=esri:102008"))
  #plot(NorthA_t)
  
  ## Create the 300km buffer zone into the ocean
  Buff <- gBuffer(NorthA_t, width = 300*1000)
  ## Remove small holes in Buff
  Buffer_t <- Polygons(Buff@polygons[[1]]@Polygons[1], ID = "Buffer")
  Buffer_t <- SpatialPolygons(list(Buffer_t), proj4string = CRS(proj4string(Buff)))
  #plot(Buffer_t)
  #plot(NorthA_t, add =T)
  
  
  ## Transform back to long lat
  Buffer <- spTransform(Buffer_t, CRS(proj4string(NorthA)))
  #plot(Buffer)
  #plot(NorthA, add = T)
  
  Buffer <- recenter(Buffer)
  NorthA <- recenter(NorthA)
  
  save(grace_raw, gia_ewh, NorthA, NorthA_t, Buffer, Buffer_t, GPS, file = datamap)
}
