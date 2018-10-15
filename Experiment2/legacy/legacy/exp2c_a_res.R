#### Assemble the results
#### Experiment 2c_a
#########################################################################################################
if(gracedata == "masscons"){
resdata <- paste0(dd,"WP1-BHM/Experiment2c/Exp2c_a.RData")
resinla <- paste0(dd, "WP1-BHM/Experiment2c/Exp2c_a.rds")
}else{
  resdata <- paste0(dd,"WP1-BHM/Experiment2c/Exp2c_a_lg.RData")
  resinla <- paste0(dd, "WP1-BHM/Experiment2c/Exp2c_a_lg.rds")
}

if(file.exists(resdata)){
  load(resdata)
}else{
  
res_inla <- readRDS(resinla)

#### Assemble the results
## Extract the predictions
INLA_massO <- res_inla$summary.random$massO
INLA_massL <- res_inla$summary.random$massL
INLA_steric <- res_inla$summary.random$steric

ssh_idx <- inla.stack.index(stkall, tag = "ssh")$data
INLA_ssh <- res_inla$summary.linear.predictor[ssh_idx,]

## mass Ocean
massO_mean <- INLA_massO$mean
massO_u <- INLA_massO$sd
proj <- inla.mesh.projector(mesh_ocean, projection = "longlat", dims = c(360,180), xlim = c(0.5,359.5), ylim = c(-89.5, 89.5))
mass_grid <- expand.grid(proj$x, proj$y)
massO_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                         mean = as.vector(inla.mesh.project(proj, as.vector(massO_mean))),
                         u = as.vector(inla.mesh.project(proj, as.vector(massO_u))))
## mass Land
massL_mean <- INLA_massL$mean
massL_u <- INLA_massL$sd
proj <- inla.mesh.projector(mesh_land, projection = "longlat", dims = c(360,180), xlim = c(0.5,359.5), ylim = c(-89.5, 89.5))
mass_grid <- expand.grid(proj$x, proj$y)
massL_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                         mean = as.vector(inla.mesh.project(proj, as.vector(massL_mean))),
                         u = as.vector(inla.mesh.project(proj, as.vector(massL_u))))
mass_pred <- na.omit(rbind(massL_pred, massO_pred))
mass_pred <- mass_pred[order(mass_pred$lat, mass_pred$lon), ]

mass_pred$mean2 <- ifelse(abs(mass_pred$mean) > 19, sign(mass_pred$mean)*20, mass_pred$mean)
mass_pred$u2 <- ifelse(mass_pred$u > 3.8, 4, mass_pred$u)
coords <- cbind(mass_pred$lon, mass_pred$lat)
mass_sp <- SpatialPoints(coords=coords, proj4string = CRS("+proj=longlat"))
mass_spg <- points2grid(mass_sp)

mass_pred$areas <- geosphere::areaPolygon(as(mass_spg, "SpatialPolygons"))/(1000^2)

## mass at GRACE resolution
if(gracedata == "masscons"){
  landmean <- drop(A_graceL %*% massL_mean)
  landerror <- drop(sqrt(A_graceL^2%*%massL_u^2))
  oceanmean <- drop(A_graceO %*% massO_mean)
  oceanerror <- drop(sqrt(A_graceO^2 %*% massO_u^2))
  grace_pred <- data.frame(id = c(goid, glid), mean = c(oceanmean, landmean), 
                         u = c(oceanerror, landerror))
  grace_pred <- grace_pred[order(grace_pred$id),]
}else{
  landmean <- drop(A_graceL %*% massL_mean)
  landsd<- drop(sqrt(A_graceL^2%*%massL_u^2))
  Land_spf$predmean <- landmean[Land_spf$id]
  Land_spf$predsd <- landsd[Land_spf$id]
  Greenland$predmean <- landmean[Greenland$id]
  Greenland$predsd <- landsd[Greenland$id]
  Antarctica$predmean <- landmean[Antarctica$id]
  Antarctica$predsd <- landsd[Antarctica$id]
  
  Ocean_spf$predmean <- drop(A_graceO %*% massO_mean)
  Ocean_spf$predsd <- drop(sqrt(A_graceO^2 %*% massO_u^2))
}


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
ssh_mean <- INLA_ssh$mean + vlm
ssh_u <- INLA_ssh$sd
ssh_pred <- data.frame(lon = steric_grid[,1], lat = steric_grid[,2],
                       mean = as.vector(inla.mesh.project(proj2, as.vector(ssh_mean))),
                       u = as.vector(inla.mesh.project(proj2, as.vector(ssh_u))))
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
## Find ids for steric 
#res_inla$misc$configs$contents
ids <- res_inla$misc$configs$contents
id1 <- ids$start[which(ids$tag == "steric")]
idn <- id1 + ids$length[which(ids$tag == "steric")]-1
stids <- id1:idn 
stericpost <- sapply(1:1000, function(x) postsamps[[x]]$latent[stids,])
bvals <- Abasin %*% stericpost
bsd <- apply(bvals, 1, sd)
gvals <- colMeans(stericpost)
Basin_pred$meanpost <- c(rowMeans(bvals), mean(gvals))

Basin_pred$upost <- c(bsd, sd(gvals))

if(gracedata == "masscons"){
ress_2c_a <- list(res_inla = res_inla, 
                  spde = list(massO = massOcean_spde, massL = massLand_spde, steric = steric_spde),  
                  mesh = list(massO = mesh_ocean, massL=mesh_land, steric = mesh_ocean),
                  pred = list(mass = mass_pred,steric = steric_pred, 
                              ssh = ssh_pred, grace_pred = grace_pred,
                              basins = Basin_pred))

save(postsamps, file = paste0(dd,"WP1-BHM/Experiment2c/Exp2c_a_postsamps.RData"))
save(ress_2c_a, file = paste0(dd,"WP1-BHM/Experiment2c/Exp2c_a.RData"))
}else{
  ress_2c_a <- list(res_inla = res_inla, 
                    spde = list(massO = massOcean_spde, massL = massLand_spde, steric = steric_spde),  
                    mesh = list(massO = mesh_ocean, massL=mesh_land, steric = mesh_ocean),
                    pred = list(mass = mass_pred,steric = steric_pred, 
                                ssh = ssh_pred, graceL = Land_spf, graceO=Ocean_spf,
                                graceGL = Greenland, graceAnt = Antarctica,
                                basins = Basin_pred))
  
  save(postsamps, file = paste0(dd,"WP1-BHM/Experiment2c/Exp2c_a_postsamps_lg2.RData"))
  save(ress_2c_a, file = paste0(dd,"WP1-BHM/Experiment2c/Exp2c_a_lg.RData"))
}




source(paste0(wd,"gmrcode/BHM_sphere/ncdfio.R"))

## save mass
df2ncdf(df = mass_pred[, c("lon", "lat", "mean", "u")], fname = paste0(dd, "WP1-BHM/Experiment2c/outputs/exp2c_a_mass_lg.nc"),  
        vars = list(c("mean", "trend", "BHM predicted mass trend", "mm/yr"), 
                    c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
        title = "BHM-predicted mass", 
        append = FALSE)

## save steric
df2ncdf(df = steric_pred[,c("lon", "lat", "mean", "u")], fname = paste0(dd, "WP1-BHM/Experiment2c/outputs/exp2c_a_steric_lg.nc"),  
        vars = list(c("mean", "trend", "BHM predicted steric trend", "mm/yr"), 
                    c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
        title = "BHM-predicted steric", 
        append = FALSE)

## save SSH
df2ncdf(df = ssh_pred[,c("lon", "lat", "mean", "u")], fname = paste0(dd, "WP1-BHM/Experiment2c/outputs/exp2c_a_ssh_lg.nc"),  
        vars = list(c("mean", "trend", "BHM predicted SSH trend", "mm/yr"), 
                    c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
        title = "BHM-predicted SSH", 
        append = FALSE)
}