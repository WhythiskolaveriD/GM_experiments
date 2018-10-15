#### Assemble the results
#### Experiment 2
#########################################################################################################

resdata <- paste0(dd_root, "WP1-BHM/Experiment2/", filename, "_res.R")
resinla <- paste0(dd_root, "WP1-BHM/Experiment2/", filename, "_inla.rds")

if(file.exists(resdata)){
  load(resdata)
}else{
  res_inla <- readRDS(resinla)
  
  #### Assemble the results
  ## Extract the predictions
  INLA_mass <- res_inla$summary.random$mass
  INLA_steric <- res_inla$summary.random$steric
  
  ## mass Ocean
  mass_mean <- INLA_mass$mean
  mass_u <- INLA_mass$sd
  proj <- inla.mesh.projector(mesh0, projection = "longlat", dims = c(360,180), xlim = c(0.5,359.5), ylim = c(-89.5, 89.5))
  mass_grid <- expand.grid(proj$x, proj$y)
  mass_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                          mean = as.vector(inla.mesh.project(proj, as.vector(mass_mean))),
                          u = as.vector(inla.mesh.project(proj, as.vector(mass_u))))
  
  mass_pred$mean2 <- ifelse(abs(mass_pred$mean) > 19, sign(mass_pred$mean)*20, mass_pred$mean)
  mass_pred$u2 <- ifelse(mass_pred$u > 3.8, 4, mass_pred$u)
  coords <- cbind(mass_pred$lon, mass_pred$lat)
  mass_sp <- SpatialPoints(coords=coords, proj4string = CRS("+proj=longlat"))
  mass_spg <- points2grid(mass_sp)
  mass_pred$areas <- geosphere::areaPolygon(as(mass_spg, "SpatialPolygons"))/(1000^2)
  
  ## steric
  steric_mean <- INLA_steric$mean
  steric_u <- INLA_steric$sd
  proj2 <- inla.mesh.projector(mesh_ocean, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-89.5, 89.5))
  steric_grid <- expand.grid(proj2$x, proj2$y)
  steric_pred <- data.frame(lon = steric_grid[,1], lat = steric_grid[,2],
                            mean = as.vector(inla.mesh.project(proj2, as.vector(steric_mean))),
                            u = as.vector(inla.mesh.project(proj2, as.vector(steric_u))))
  steric_pred$mean2 <- ifelse(abs(steric_pred$mean) > 19, sign(steric_pred$mean)*20, steric_pred$mean)
  steric_pred$u2 <- ifelse(steric_pred$u > 3.8, 4, steric_pred$u)
  
  
  ## ssh
  ssh_mean <- steric_pred$mean + mass_pred$mean
  ssh_u <- sqrt(steric_pred$u^2 + mass_pred$u^2)
  ssh_pred <- data.frame(lon = steric_grid[,1], lat = steric_grid[,2],
                         mean = ssh_mean, u = ssh_u)
  ssh_pred$mean2 <- ifelse(abs(ssh_pred$mean) > 19, sign(ssh_pred$mean)*20, ssh_pred$mean)
  ssh_pred$u2 <- ifelse(ssh_pred$u > 3.8, 4, ssh_pred$u)
  
  ## Basins
  ## Use the MAP estimator of the parameters to get the covariance matrix of steric
  theta_post <- summary(res_inla)$hyperpar
  Q <- inla.spde.precision(steric_spde, theta = theta_post$mode[5:6])
  V <- cov2cor(inla.qinv(Q))
  Basin_mean <- drop(Abasin %*% steric_mean)
  Global_mean <- mean(na.omit(steric_mean))
  Basin_sd <- drop(sqrt(diag(Abasin %*% (Diagonal(x = steric_u) %*% V %*% Diagonal(x = steric_u)) %*% t(Abasin))))
  Global_sd <- sqrt(sum((Diagonal(x = steric_u) %*% V %*% Diagonal(x = steric_u)))/nrow(V)^2)
  Basin_sd2 <- drop(sqrt(diag(Abasin %*% (Diagonal(x = steric_u)^2) %*% t(Abasin))))
  Global_sd2 <- sqrt(sum(steric_u^2)/nrow(V)^2)
  Basin_pred <- data.frame(mean = c(Basin_mean, Global_mean), u = c(Basin_sd, Global_sd),
                           u2 = c(Basin_sd2, Global_sd2))
  
  ## Samples form posteriors 
  postsamps <- inla.posterior.sample(n = 1000, result=res_inla, use.improved.mean = TRUE, add.names = FALSE, seed = 12)
  ## Find ids for the latent process
  #res_inla$misc$configs$contents
  ids <- res_inla$misc$configs$contents
  ids$end <- ids$start + ids$length -1
  massids <- ids$start[which(ids$tag == "mass")] : ids$end[which(ids$tag == "mass")]
  stericids <- ids$start[which(ids$tag == "steric")] : ids$end[which(ids$tag == "steric")]
  masspost <- sapply(1:1000, function(x) postsamps[[x]]$latent[massids,])
  stericpost <- sapply(1:1000, function(x) postsamps[[x]]$latent[stericids,])
  ## calculate posterior for ssh
  
  
  
  bvals <- Abasin %*% stericpost
  bsd <- apply(bvals, 1, sd)
  gvals <- colMeans(stericpost)
  Basin_pred$meanpost <- c(rowMeans(bvals), mean(gvals))
  
  Basin_pred$upost <- c(bsd, sd(gvals))
  
  ress_2 <- list(res_inla = res_inla, 
                 spde = list(mass = mass_spde, steric = steric_spde),  
                 mesh = list(mass = mesh0, steric = mesh_ocean),
                 pred = list(mass = mass_pred,steric = steric_pred, 
                             ssh = ssh_pred, basins = Basin_pred))
  
  save(postsamps, file = paste0(dd_root,"WP1-BHM/Experiment2/", filename, "_postsamps.RData"))
  save(ress_2, file = paste0(dd_root,"WP1-BHM/Experiment2/", filename, "_res.RData"))
  
  
  ########################################################################################
  ### Save ncdf files
  
  ## save mass
  df2ncdf(df = mass_pred[, c("lon", "lat", "mean", "u")], 
          fname = paste0(dd_root, "WP1-BHM/Experiment2/outputs/", filename, "_mass.nc"),  
          vars = list(c("mean", "trend", "BHM predicted mass trend", "mm/yr"), 
                      c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
          title = "BHM-predicted mass", 
          append = FALSE)
  
  ## save steric
  df2ncdf(df = steric_pred[,c("lon", "lat", "mean", "u")],
          fname = paste0(dd_root, "WP1-BHM/Experiment2/outputs/", filename, "_steric.nc"),  
          vars = list(c("mean", "trend", "BHM predicted steric trend", "mm/yr"), 
                      c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
          title = "BHM-predicted steric", 
          append = FALSE)
  
  ## save SSH
  df2ncdf(df = ssh_pred[,c("lon", "lat", "mean", "u")], 
          fname = paste0(dd_root, "WP1-BHM/Experiment2/outputs/", filename, "_ssh.nc"),  
          vars = list(c("mean", "trend", "BHM predicted SSH trend", "mm/yr"), 
                      c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
          title = "BHM-predicted SSH", 
          append = FALSE)
  
}

## Add group read and write permissions to all the files recursively to the experiment 2 folder
if(Sys.info()["nodename"] != "IT034844"){
  system(paste0("chmod -R g+rw ", dd_root, "WP1-BHM/Experiment2/"))
}

## You can also add the permissions to the files by using
## Sys.chmod(path2file, mode = "0770", use_umask = FALSE)
## Note that use_umask for the rdsf is more restricted than the working dir, so need to set use_umask to FALSE
