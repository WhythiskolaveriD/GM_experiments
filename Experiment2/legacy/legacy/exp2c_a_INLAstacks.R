#### Build INLA spdes and stacks
#### Experiment 2c_a
spdes <- paste0(dd,"WP1-BHM/Experiment2c/data/exp2c_spdes.rda")
if(gracedata == "masscons"){
  stks <- paste0(dd,"WP1-BHM/Experiment2c/data/exp2c_staks.rda")
}else{
  stks <- paste0(dd,"WP1-BHM/Experiment2c/data/exp2c_staks_lg.rda")
}


#### Build SPDES
#########################################################################################################
if(file.exists(spdes)){
  load(spdes)
}else{
  massOcean_spde <- inla.spde2.matern(mesh_ocean, B.tau = matrix(c(ltaug1, -1, 1),1,3), B.kappa = matrix(c(lkappag1, 0, -1), 1,3),
                                      theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/thetag1_ss), sqrt(1/thetag1_rs)))
  massLand_spde <- inla.spde2.matern(mesh_land, B.tau = matrix(c(ltaug2, -1, 1),1,3), B.kappa = matrix(c(lkappag2, 0, -1), 1,3),
                                     theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/thetag2_ss), sqrt(1/thetag2_rs)))
  steric_spde <- inla.spde2.matern(mesh_ocean, B.tau = matrix(c(ltau2, -1, 1),1,3), B.kappa = matrix(c(lkappa2, 0, -1), 1,3),
                                   theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta2_ss), sqrt(1/theta2_rs)))
  save(massOcean_spde, massLand_spde, steric_spde, file = spdes)
}


#### Build the stacks
#########################################################################################################
if(file.exists(stks)){
  load(stks)
}else{
  ## The altimetry data
  yaltv <- alt_data2$trend_vlm
  ## Link Altimetry
  A_altv1 <- inla.spde.make.A(mesh = mesh_ocean, loc = alt_loc)
  A_altv2 <- inla.spde.make.A(mesh = mesh_ocean, loc = alt_loc)
  stk_alt <- inla.stack(data = list(y=yaltv), A = list(A_altv1,  A_altv2),
                        effects = list(list(massO = 1:massOcean_spde$n.spde),
                                       list(steric = 1:steric_spde$n.spde)), tag = "alt")
  ## Create the prediction stacks
  A_ssh <- Diagonal(n= steric_spde$n.spde, x =1)
  stk_predSSH <- inla.stack(data = list(y=NA), A = list(A_ssh, A_ssh),
                            effects = list(list(massO = 1:massOcean_spde$n.spde),
                                           list(steric = 1:steric_spde$n.spde)), tag = "ssh")
  if(gracedata == "masscons"){
    ## The grace data
    ygrace_Ocean <- grace_sp$trendgia[goid]
    ygrace_Land <- grace_sp$trendgia[-goid]
    ## Link GRACE
    grace_locO <- do.call(cbind, Lll2xyz(lat = grace_sp@coords[goid,2], lon = grace_sp@coords[goid,1]))
    grace_locL <- do.call(cbind, Lll2xyz(lat = grace_sp@coords[-goid,2], lon = grace_sp@coords[-goid,1]))
    A_graceO <- inla.spde.make.A(mesh = mesh_ocean, loc = grace_locO)
    A_graceL <- inla.spde.make.A(mesh = mesh_land, loc = grace_locL)
    ## Create the data stacks
    stk_graceO <- inla.stack(data = list(y=ygrace_Ocean), A = list(A_graceO),
                             effects = list(massO = 1:massOcean_spde$n.spde), tag = "graceO")
    
    stk_graceL <- inla.stack(data = list(y=ygrace_Land), A = list(A_graceL),
                             effects = list(massL = 1:massLand_spde$n.spde), tag = "graceL")
    ## prediction of the the grace
    ## Dont! use projection afterwards!
    stkall <- inla.stack(stk_graceO, stk_graceL, stk_alt, stk_predSSH)
    prec_scale <- c(1/grace_sp$std[goid]^2, 1/grace_sp$std[-goid]^2, 1/alt_data2$err_ssh^2, 
                    rep(100, steric_spde$n.spde))
  }else{
    ## First create the block ID for each of the mesh points
    ## Land
    Landll<- Lxyz2ll(list(x=mesh_land$loc[,1], y = mesh_land$loc[,2], z = mesh_land$loc[,3]))
    land_sp <- SpatialPoints(coords = cbind(Landll$lon, Landll$lat), proj4string = CRS(proj4string(Land_spf)))
    land_sp2 <- spTransform(land_sp, proj4string(Greenland))
    land_sp3 <- spTransform(land_sp, proj4string(Antarctica))
    ## Add an id column to the catchment region
    Land_spf$id <- 1:nrow(Land_spf@data)
    Greenland$id <- (1+nrow(Land_spf@data)):(nrow(Land_spf@data)+nrow(Greenland@data))
    Antarctica$id <- (1+Greenland$id[length(Greenland$id)]):(Greenland$id[length(Greenland$id)] + nrow(Antarctica@data)) 
    land_id <- over(land_sp, Land_spf)$id
    green_id <- over(land_sp2, Greenland)$id 
    Ant_id <- over(land_sp3, Antarctica)$id
    all_id <- land_id
    all_id[which(!is.na(green_id))] <- green_id[which(!is.na(green_id))]
    all_id[which(!is.na(Ant_id))] <- Ant_id[which(!is.na(Ant_id))]
    ## Find points inside the catchamet region -- id is not NA
    landloc <- mesh_land$loc[which(!is.na(all_id)),]
    landblock <- as.numeric(na.omit(all_id))
    ## Build the link matrix
    A_graceL <- inla.spde.make.A(mesh = mesh_land, loc = landloc, block = landblock, block.rescale = "count" )
    ## Build the data stack
    stk_graceL <- inla.stack(data = list(y=c(Land_spf$trend, Greenland$trend, Antarctica$trend)), A = list(A_graceL),
                             effects = list(list(massL = 1:massLand_spde$n.spde)), tag = "graceL")
    
    ## Ocean
    Oceanll<- Lxyz2ll(list(x=mesh_ocean$loc[,1], y = mesh_ocean$loc[,2], z = mesh_ocean$loc[,3]))
    Oceanll$lon <- ifelse(Oceanll$lon < 0, Oceanll$lon + 360, Oceanll$lon)
    ocean_sp <- SpatialPoints(coords = cbind(Oceanll$lon, Oceanll$lat), proj4string = CRS(proj4string(Ocean_spf)))
    ## Add an id column to the ocean polygons
    Ocean_spf$id <- 1:nrow(Ocean_spf@data)
    ocean_id <- over(ocean_sp, Ocean_spf)$id
    oceanloc <- mesh_ocean$loc[which(!is.na(ocean_id)),]
    oceanblock <- na.omit(ocean_id)
    ## Build the link matrix
    A_graceO <- inla.spde.make.A(mesh = mesh_ocean, loc = oceanloc, block = oceanblock, block.rescale = "count" )
    ## Build the data stack
    stk_graceO <- inla.stack(data = list(y=Ocean_spf$trend), A = list(A_graceO),
                             effects = list(list(massO = 1:massOcean_spde$n.spde)), tag = "graceO")
    
    stkall <- inla.stack(stk_graceO, stk_graceL, stk_alt, stk_predSSH)
    
    prec_scale1 <- c(1/Ocean_spf$std^2, 1/Land_spf$std^2, 1/Greenland$std^2, 1/Antarctica$std^2,
                    1/alt_data2$err_ssh^2, rep(100, steric_spde$n.spde))
    prec_scale2 <- c(1/Ocean_spf$std2^2, 1/Land_spf$std2^2, 1/Greenland$std2^2, 1/Antarctica$std2^2,
                     1/alt_data2$err_ssh^2, rep(100, steric_spde$n.spde))
  }

  ## Predict Basin level steric
  ## Load the Basin polygons
  load(paste0(dd, "WP1-BHM/maps/Ocean/basins.rda"))
  ## Number of Basins
  nrow(allbasinsdf)
  ## convert the mesh vertice location to longlat
  vll <- Lxyz2ll(list(x=mesh_ocean$loc[,1], y = mesh_ocean$loc[,2], z = mesh_ocean$loc[,3]))
  vll$lon <- ifelse(vll$lon < 0, vll$lon + 360, vll$lon)
  v_sp <- SpatialPoints(coords = cbind(vll$lon, vll$lat))
  proj4string(v_sp) <- proj4string(allbasinsdf) <- CRS("+proj=longlat")
  BasinId <- over(v_sp, allbasinsdf)
  Bnames <- levels(BasinId$names)
  Basin_loc <- mesh_ocean$loc[which(!is.na(BasinId$names)),]
  Basin_block <- BasinId$names[which(!is.na(BasinId$names))]
  Basin_bkid <- as.numeric(Basin_block)
  ## Build the matrix that link the process to the Basin
  Abasin <- inla.spde.make.A(mesh=mesh_ocean, loc = Basin_loc, block = Basin_bkid, block.rescale = "count")
  
  ## save the objects
  if(gracedata == "masscons"){
    save(stkall, prec_scale, A_graceL, A_graceO, Abasin, allbasinsdf,  file = stks)
  }else{
    save(stkall, prec_scale1, prec_scale2, A_graceL, A_graceO, Abasin, allbasinsdf,  file = stks)
  }

  
}


