###### This script generate the mesh and processing data for Experiment 3
####################################################################################################
meshdata <- paste0(dd_save, expname, "_meshdata.rda")

if(file.exists(meshdata)){
  load(meshdata)
}else{
  #### Generate mesh
  ####################################################################################################
  ## sample points over the region (1) evenly distributed (2) dense at the boundary
  ## Try one degree resolution first -- about 2500 cells 
  Regpoints_t <- spsample(Buffer_t, n = 3000, type = "regular")
  #plot(Buffer_t)
  #plot(NorthA_t, add = T)
  #points(Regpoints_t@coords, pch = ".")
  
  ## Transform to long lat
  Regpoints <- spTransform(Regpoints_t, CRS(proj4string(NorthA)))
  Regpoints@coords[,1] <- ifelse(Regpoints@coords[,1] < 0, Regpoints@coords[,1] + 360, Regpoints@coords[,1])
  #plot(Buffer)
  #plot(NorthA, add =T)
  #points(Regpoints@coords, pch = ".")
  #points(GPS@coords, col = 2, pch = "*")
  
  ## Add more points along boundary 
  ## create buffer zones along the boundary
  buff_t1 <- gBuffer(NorthA_t, width = -100*1000)
  buff_t2 <- gBuffer(NorthA_t, width =  100*1000)
  buff_tb <- gDifference(buff_t2, buff_t1)
  bps_t <- spsample(buff_tb, n = 3000, type = "regular")
  
  bps <- spTransform(bps_t, CRS(proj4string(NorthA)))
  bps@coords[,1] <- ifelse(bps@coords[,1] < 0, bps@coords[,1] + 360, bps@coords[,1])
  
  #plot(Buffer_t)
  #plot(NorthA_t, add = T)
  #points(Regpoints_t@coords, pch = ".")
  #points(bps_t@coords, pch = ".", col = 2, cex = 2)
  
  #plot(Buffer)
  #plot(NorthA, add = T)
  #points(Regpoints@coords, pch = ".")
  #points(bps@coords, pch = ".", col = 2, cex = 2)
  
  ## conbine all the points and transform to xyz
  mesh_ll <- rbind(Regpoints@coords, bps@coords)
  mesh_xyz <- do.call(cbind, GEOmap::Lll2xyz(lat = mesh_ll[,2], lon = mesh_ll[,1]))
  mesh0 <- inla.mesh.2d(loc = mesh_xyz, cutoff = 0.005, max.edge = 0.05)
  #plot(mesh0)
  
  ## Extract the outer boundary of the mesh
  mesh_bnd <- inla.mesh.boundary(mesh0)
  mesh_boundll <- GEOmap::Lxyz2ll(list(x = mesh_bnd[[1]]$loc[,1], y = mesh_bnd[[1]]$loc[,2], z = mesh_bnd[[1]]$loc[,3]))
  mesh_boundll$lon <- ifelse(mesh_boundll$lon < 0, mesh_boundll$lon + 360, mesh_boundll$lon)
  mesh_boundll <- cbind(mesh_boundll$lon, mesh_boundll$lat)[mesh_bnd[[1]]$idx[,1],]
  mb_poly <- Polygon(coords = mesh_boundll)
  mb_poly <- Polygons(list(mb_poly), ID = "mb")
  mb_poly <- SpatialPolygons(list(mb_poly), proj4string = CRS(proj4string(NorthA)))
  
  #### Extract raw data in this region
  ####################################################################################################
  ## GPS data -- all within this region
  #plot(mb_poly)
  #points(GPS@coords)
  
  ## GIA data
  proj4string(gia_ewh) <- CRS(proj4string(mb_poly))
  gia_mb <- which(over(gia_ewh, mb_poly)==TRUE)
  gia_northA <- gia_ewh[gia_mb, ]
 
  ## GRACE data
  proj4string(grace_raw) <- CRS(proj4string(mb_poly))
  grace_mb <- which(over(grace_raw, mb_poly)==TRUE)
  grace_northA <- grace_raw[grace_mb,]
  
  ## Create the GIA prior at the mesh vertices
  ## Find which grid cell of gia the mesh vertice fall
  gia_NAg <- gia_ewh
  gridded(gia_NAg) <- TRUE
  meshll <- GEOmap::Lxyz2ll(list(x = mesh0$loc[,1], y = mesh0$loc[,2], z = mesh0$loc[,3]))
  meshll$lon <- ifelse(meshll$lon < 0, meshll$lon + 360, meshll$lon)
  meshspp <- SpatialPoints(coords = cbind(meshll$lon, meshll$lat), proj4string = CRS(proj4string(gia_northA)))
  meshgia <- over(meshspp, gia_NAg)$trend
  
  ## detrend the GPS observation
  proj4string(GPS) <- proj4string(gia_NAg)
  gpsgia <- GPS$trend - over(GPS, gia_NAg)$trend
  GPS$detrend <- gpsgia
  
  ## detrend the GRACE observation
  gracegia <- grace_northA$trend - over(grace_northA, gia_NAg)$trend
  grace_northA$detrend <- gracegia
  
  save(mesh0, mb_poly, meshgia, gia_northA, grace_northA, GPS, file = meshdata)
}
