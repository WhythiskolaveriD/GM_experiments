#### Build INLA spdes and stacks
stks <- paste0(dd_root,"WP1-BHM/Experiment2/data/Exp2_GSFCregion_stks.rda")

## Note that this script has the MAJOR DIFFERENCE from the GSFCgrid script!!!
## Here the process points are linked to regional observations!!

if(file.exists(stks)){
  load(stks)
}else{
  #### Build the stacks
  #########################################################################################################
  ## The altimetry data
  yaltv <- alt_data2$trend
  ## Link Altimetry
  A_altv2Mass <- inla.spde.make.A(mesh = mesh0, loc = alt_loc)
  A_altv2Steric <- inla.spde.make.A(mesh = mesh_ocean, loc = alt_loc)
  stk_alt <- inla.stack(data = list(y=yaltv), A = list(A_altv2Mass,  A_altv2Steric),
                        effects = list(list(mass = 1:mass_spde$n.spde),
                                       list(steric = 1:steric_spde$n.spde)), tag = "alt")
  
  ## The grace data
  ## First create the block ID for each of the mesh points 
  ## The lonlat of each mesh points are already saved in the object "mesh_locll" in the previous script
  mesh_sp <- SpatialPoints(coords = cbind(mesh_locll[,1], mesh_locll[,2]), proj4string = CRS(proj4string(Land_spf)))
  
  ## Create ID for Land
  mesh_sp2 <- spTransform(mesh_sp, proj4string(Greenland))
  mesh_sp3 <- spTransform(mesh_sp, proj4string(Antarctica))
  
  ## Add an id column to the catchment region
  Land_spf$id <- 1:nrow(Land_spf@data)
  Greenland$id <- (1+nrow(Land_spf@data)):(nrow(Land_spf@data)+nrow(Greenland@data))
  Antarctica$id <- (1+Greenland$id[length(Greenland$id)]):(Greenland$id[length(Greenland$id)] + nrow(Antarctica@data)) 
  land_id <- over(mesh_sp, Land_spf)$id
  green_id <- over(mesh_sp2, Greenland)$id 
  Ant_id <- over(mesh_sp3, Antarctica)$id
  
  land_id[which(!is.na(green_id))] <- green_id[which(!is.na(green_id))]
  land_id[which(!is.na(Ant_id))] <- Ant_id[which(!is.na(Ant_id))]
  landloc <- mesh0$loc[which(!is.na(land_id)),]
  landid <- as.numeric(na.omit(land_id))
  
  ## Create ID for Ocean
  mesh_sp4 <- mesh_sp
  mesh_sp4@coords[,1] <-  ifelse(mesh_sp4@coords[,1] < 0, mesh_sp4@coords[,1] + 360, mesh_sp4@coords[,1])
  proj4string(mesh_sp4) <- CRS(proj4string(Ocean_spf))
  
  ## Add an id column to the ocean polygons
  Ocean_spf$id <- (max(landid)+1):(max(landid)+nrow(Ocean_spf@data))
  ocean_id <- over(mesh_sp4, Ocean_spf)$id
  oceanloc <- mesh0$loc[which(!is.na(ocean_id)),]
  oceanid <- na.omit(ocean_id)
  
  all_loc <- rbind(landloc, oceanloc)
  allblock <- c(landid, oceanid)
  
  A_grace <- inla.spde.make.A(mesh = mesh0, loc = all_loc, block = allblock, block.rescale = "count" )
  
  ## Build the data stack
  stk_grace <- inla.stack(data = list(y=c(Land_spf$trend, Greenland$trend, Antarctica$trend, Ocean_spf$trend)), A = list(A_grace),
                           effects = list(list(massL = 1:mass_spde$n.spde)), tag = "grace")
  
  
  stkall <- inla.stack(stk_grace, stk_alt)
  
  
  prec_scale <- c(1/Land_spf$std^2, 1/Greenland$std^2, 1/Antarctica$std^2, 1/Ocean_spf$std^2, 1/alt_data2$err_ssh^2) # inflated error
  prec_scale1 <- c(1/Land_spf$std1^2, 1/Greenland$std1^2, 1/Antarctica$std1^2, 
                   1/Ocean_spf$std1^2, 1/alt_data2$err_ssh^2) # inflated error + gia error
  
  ## Predict Basin level steric
  ## Load the Basin polygons
  load(paste0(dd_root, "WP1-BHM/maps/Ocean/basins.rda"))
  ## Number of Basins
  nrow(allbasinsdf)
  ## convert the mesh vertice location to longlat
  vll <- Lxyz2ll(list(x=mesh0$loc[,1], y = mesh0$loc[,2], z = mesh0$loc[,3]))
  vll$lon <- ifelse(vll$lon < 0, vll$lon + 360, vll$lon)
  v_sp <- SpatialPoints(coords = cbind(vll$lon, vll$lat))
  proj4string(v_sp) <- proj4string(allbasinsdf) <- CRS("+proj=longlat")
  BasinId <- over(v_sp, allbasinsdf)
  Bnames <- levels(BasinId$names)
  Basin_loc <- mesh0$loc[which(!is.na(BasinId$names)),]
  Basin_block <- BasinId$names[which(!is.na(BasinId$names))]
  Basin_bkid <- as.numeric(Basin_block)
  ## Build the matrix that link the process to the Basin
  Abasin <- inla.spde.make.A(mesh=mesh_ocean, loc = Basin_loc, block = Basin_bkid, block.rescale = "count")
  
  save(stkall, prec_scale, prec_scale1, A_grace, Abasin, allbasinsdf,  file = stks)
  
}

