#### Build INLA spdes and stacks
stks <- paste0(dd_root,"WP1-BHM/Experiment2/data/Exp2_GSFCgrid_stks.rda")

## The process point is linked to the the grid observations.
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
  ygrace <- grace_sp$trendgia
  
  ## Link GRACE
  grace_loc <- do.call(cbind, Lll2xyz(lat = grace_sp@coords[,2], lon = grace_sp@coords[,1]))
  A_grace <- inla.spde.make.A(mesh = mesh0, loc = grace_loc)
  
  ## Create the data stacks
  stk_grace <- inla.stack(data = list(y=ygrace), A = list(A_grace),
                          effects = list(mass = 1:mass_spde$n.spde), tag = "grace")
  
  stkall <- inla.stack(stk_grace, stk_alt)
  prec_scale <- c(1/grace_sp$std^2,  1/alt_data2$err_ssh^2) # raw error
  prec_scale1 <- c(1/grace_sp$std1^2,  1/alt_data2$err_ssh^2) # inflated error
  prec_scale2 <- c(1/grace_sp$std2^2,  1/alt_data2$err_ssh^2) # inflated error + gia error
  
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
  
  save(stkall, prec_scale, prec_scale1, prec_scale2, A_grace, Abasin, allbasinsdf,  file = stks)
  
}

