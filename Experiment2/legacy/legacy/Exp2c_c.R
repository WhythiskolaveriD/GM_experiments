## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

if(Sys.info()["nodename"] == "IT034844"){
  ## BS
  wd <- "C:/ZSwork/"
  dd <- "Z:/"
}else if(Sys.info()["nodename"] == "IT064613"){
  ## Maths
  wd <-"/home/zs16444/coding/"
  dd <-"/home/zs16444/globalmass/"
}else{
  ## server  
  wd <- "~/"
  dd <- "/./projects/GlobalMass/"
}

library(rgdal); library(sp);library(GEOmap)
library(INLA)
library(ncdf4)
source(paste0(wd,"gmrcode/BHM_sphere/functions.R"))
source(paste0(wd, "gmrcode/BHM_sphere/partition_fun.R"))
source(paste0(wd, "gmrcode/BHM_sphere/pp_functions.R"))

## ----data----------------------------------------------------------------
source(paste0(wd,"experiments/Experiments2/Doc/c/data_Ex2c_c.R"))

## ----ssh_mesh------------------------------------------------------------
## Genereate Fibonacci points on the sphere
fibo_points <- fiboSphere(N = 3e4, LL = FALSE)
mesh0 <- inla.mesh.2d(loc = fibo_points, cutoff = 0.01, max.edge = 1)
summary(mesh0)

## ----ssh_mesh2-----------------------------------------------------------
mesh0 <- dt.mesh.addon.posTri(mesh = mesh0, globe = TRUE) 
Tlonlat <- Lxyz2ll(list(x = mesh0$posTri[,1], y = mesh0$posTri[,2], z = mesh0$posTri[,3]))
mesh0$Trill <- cbind(lon = Tlonlat$lon, lat =Tlonlat$lat)
T_sp <- SpatialPoints(coords=mesh0$Trill, proj4string = CRS(proj4string(Land)))
## Devide these points into three groups ( > +- 85)
Npole_id <- which(T_sp@coords[,2] > 85)
Spole_id <- which(T_sp@coords[,2] < -85)
mid_id <- (1:nrow(T_sp@coords))[-c(Npole_id, Spole_id)]
Spole_sp <- T_sp[Spole_id]
mid_sp <- T_sp[mid_id]
## Nothing to be done on the north pole region points since they will all be classified in the Ocean.

## Transform the middle part and South pole region into suitable projection system
mid_sp <- spTransform(mid_sp, CRS("+init=epsg:3857"))
Spole_sp <- spTransform(Spole_sp, CRS("+init=epsg:5481"))

## Create the Land label
TinLand1 <- mid_id[unlist(over(Land1, mid_sp, returnList=T))]
TinLand2 <- Spole_id[unlist(over(Land2, Spole_sp, returnList=T))]
TinLand <- c(TinLand1, TinLand2)
plot(T_sp[TinLand], pch = ".")

## Create the independent zone label
TinID1 <- mid_id[unlist(over(Indzone1, mid_sp, returnList=T))]
TinID2 <- mid_id[unlist(over(Indzone2, Spole_sp, returnList=T))]
TinID <- c(TinID1, TinID2)
plot(T_sp[TinID], pch = ".")

## Finally the Ocean label
TinOcean <- (1:nrow(T_sp@coords))[-c(TinLand, TinID)]
plot(T_sp[TinOcean], pch = ".", col = 2, add = T)

## check plot
plot(T_sp, pch = ".")
plot(T_sp[TinOcean], pch = ".", col = 2, add = T)
plot(T_sp[TinLand], pch = ".", col = 3, add = T)
plot(T_sp[TinID], pch = ".", col = 4, add = T)

Omega = dt.Omega(list(TinLand, TinID, TinOcean), mesh0)
plot(mesh0, t.sub = Omega[[1]], rgl = TRUE, col = "yellow")
plot(mesh0, t.sub = Omega[[2]], rgl = TRUE, col = "grey", add = TRUE)
plot(mesh0, t.sub = Omega[[3]], rgl = TRUE, col = "blue", add = TRUE)

#play3d(spin3d(axis = c(1,1,1)), duration = 8)

movie3d( spin3d(axis = c(1,1,1)), duration = 8 )


par3dinterp(time = 0:1, userMatrix = list(M,rotate3d(M, pi/2, 0.5, 0, 0))

## ----spde----------------------------------------------------------------
## mass
massQ <- dt.create.Q(mesh0, Omega,  fixed.ranges = c(NA, 0.01, NA))
## Transform the parameters for the SPDE_GMRF approximation
#prior <- list(sigma = rbind(tsigmag1, tsigmag2), rho = rbind(trhog1, trhog2))
prior <- list(sigma = tsigmag1, range = rbind(trhog1, trhog2))
log.prior <- dt.create.prior.log.norm(prior.param = prior) 
mass_spde <- pp.inla.model(Q = massQ, log.prior=log.prior)

## steric
steric_spde <- inla.spde2.matern(mesh_ocean, B.tau = matrix(c(ltau2, -1, 1),1,3), B.kappa = matrix(c(lkappa2, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta2_ss), sqrt(1/theta2_rs)))

## ----link_data-----------------------------------------------------------
## The data
ygrace <- grace_df$trendgia
yaltv<- alt_data2$trend_vlm

## The A matrice for linking observation and processes

## Link GRACE to mass
A_grace <- inla.spde.make.A(mesh = mesh0, loc = grace_loc)

## Link Altimetry and grace to steric and mass
A_altv_s <- inla.spde.make.A(mesh = mesh_ocean, loc = alt_loc)
A_altv_m <- inla.spde.make.A(mesh = mesh0, loc = alt_loc)

## The A matrices for predict ssh
A_ssh_s <- A_ssh_m <- Diagonal(n= steric_spde$n.spde, x =1)
A_ssh_m <- inla.spde.make.A(mesh = mesh0, loc = mesh_ocean$loc)

## The errors
prec_scale <- c(1/grace_df$std^2, 1/alt_data2$err_ssh^2, rep(1, nrow(A_ssh_s)))

## Create the mass and steric stacks
stkmass <- inla.stack(data = list(y=ygrace), A = list(A_grace),
                     effects = list(mass = 1:(mesh0$n)), tag = "mass")

stksteric <- inla.stack(data = list(y=yaltv), A = list(A_altv_m, A_altv_s),
                     effects = list(list(mass = 1:mesh0$n),
                                    list(steric = 1:steric_spde$n.spde)), tag = "steric")

## Create the stack for predicting SSH
stkssh <- inla.stack(data = list(y=NA), A = list(A_ssh_m, A_ssh_s),
                     effects = list(list(mass = 1:mesh0$n), 
                                    list(steric = 1:steric_spde$n.spde)), 
                     remove.unused = FALSE, tag = "ssh")
## join the two stacks
stkall <- inla.stack(stkmass,  stksteric, stkssh, remove.unused = FALSE)

## ----inla_run, include = TRUE, eval = FALSE------------------------------
## ## Fix altimetry errors as they are known
## hyper <- list(prec = list(fixed = TRUE, initial = 0))
## 
## ## The formular -- we add the constraint that mass change sum to zero
## ## constraint 1 -- vertices sum to zero
## A1 <- matrix(1, nrow = 1, ncol = mesh0$n)
## ## Constraint 2 -- grace loc sum to zero
## A2 <- matrix(colSums(A_grace), nrow = 1)
## ## Constraint 3 -- longlat grid sum to zero
## gridll <- ice6g_vlm@coords
## gridxyz <- do.call(cbind, Lll2xyz(lat = gridll[,2], lon = gridll[,1]))
## A_grid <- inla.spde.make.A(mesh = mesh0, loc = gridxyz)
## weights <- ice6g_vlm$areas/sum(ice6g_vlm$areas)
## A3 <- matrix(weights, nrow = 1) %*% A_grid
## A <- as.matrix(rbind(A1, A2, A3))
## ## The formular
## formula = y ~ -1 + f(mass, model = mass_spde,
##                      extraconstr = list(A = A, e = c(0, 0, 0))) + f(steric, model = steric_spde)
## 
## ## Run INLA
## 
## res_inla <- inla(formula, data = inla.stack.data(stkall), family = "gaussian",
##                  scale =prec_scale, control.family = list(hyper = hyper),
##                  control.predictor=list(A=inla.stack.A(stkall), compute =TRUE), verbose = TRUE)

## ----inla_res, include = TRUE, eval = FALSE, echo = FALSE----------------
## ## idices of the predicted processes
## mass_idx <- inla.stack.index(stkall, tag = "mass")$effects
## mass_idxg <- inla.stack.index(stkall, tag = "mass")$data
## steric_idx <- inla.stack.index(stkall, tag = "steric")$effects[-c(1:mesh0$n)]
## ssh_idx <- inla.stack.index(stkall, tag = "ssh")$data
## 
## 
## ## Extract the predictions
## INLA_mass <- res_inla$summary.random$mass
## INLA_steric <- res_inla$summary.random$steric
## INLA_ssh <- res_inla$summary.linear.predictor[ssh_idx,]
## INLA_grace <- res_inla$summary.linear.predictor[mass_idxg,]
## 
## ## mass
## ## mass Ocean
## mass_mean <- INLA_mass$mean
## mass_u <- INLA_mass$sd
## proj <- inla.mesh.projector(mesh0, projection = "longlat", dims = c(360,180), xlim = c(0.5,359.5), ylim = c(-89.5, 89.5))
## mass_grid <- expand.grid(proj$x, proj$y)
## mass_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
##                        mean = as.vector(inla.mesh.project(proj, as.vector(mass_mean))),
##                        u = as.vector(inla.mesh.project(proj, as.vector(mass_u))))
## 
## coords <- cbind(mass_pred$lon, mass_pred$lat)
## mass_sp <- SpatialPoints(coords=coords, proj4string = CRS("+proj=longlat"))
## mass_spg <- points2grid(mass_sp)
## grid_areas <- geosphere::areaPolygon(as(mass_spg, "SpatialPolygons"))/(1000^2)
## mass_pred$area <- grid_areas
## 
## 
## ## mass at GRACE resolution
## grace_pred <- data.frame(mean = INLA_grace$mean, u = INLA_grace$sd)
## 
## ## steric
## steric_mean <- INLA_steric$mean
## steric_u <- INLA_steric$sd
## proj2 <- inla.mesh.projector(mesh_ocean, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-89.5, 89.5))
## steric_grid <- expand.grid(proj2$x, proj2$y)
## steric_pred <- data.frame(lon = steric_grid[,1], lat = steric_grid[,2],
##                        mean = as.vector(inla.mesh.project(proj2, as.vector(steric_mean))),
##                        u = as.vector(inla.mesh.project(proj2, as.vector(steric_u))))
## 
## ## ssh
## ssh_mean <- INLA_ssh$mean + vlm
## ssh_u <- INLA_ssh$sd
## ssh_pred <- data.frame(lon = steric_grid[,1], lat = steric_grid[,2],
##                        mean = as.vector(inla.mesh.project(proj2, as.vector(ssh_mean))),
##                        u = as.vector(inla.mesh.project(proj2, as.vector(ssh_u))))
## 
## ress_2c_c <- list(res_inla = res_inla, st = stkall,
##                 spde = list(mass = mass_spde, steric = steric_spde),
##                 mesh = list(mass = mesh0, steric = mesh_ocean),
##                 pred = list(mass = mass_pred, steric = steric_pred, ssh = ssh_pred, grace_pred = grace_pred))
## 
## save(ress_2c_c, file = "/./projects/GlobalMass/WP1-BHM/Experiment2c/Exp2c_c.RData")
## 

## ----hyper, include=TRUE, echo = FALSE-----------------------------------
load(paste0(dd, "WP1-BHM/Experiment2c/Exp2c_c.RData")) 

pars_mass <- marginal_par(res = ress_2c_c, process = "mass", SPDE2 = FALSE, 
                          theta.names = c("Sigma","rho_L","rho_O"), plot = TRUE) 

pars_steric <- marginal_par(res = ress_2c_c, process = "steric", SPDE2 = TRUE, 
                            theta.names = c("rho_S","sigma_S"), plot = TRUE) 

## The posterior modes
print(paste("The estimated correlation length for mass in Ocean is:", pars_mass$thetam[3]*6371, "km", sep = " "))
#print(paste("The estimated marginal variance for mass in Ocean is:", sqrt(pars_mass$thetam[3]), "mm/yr", sep = " "))

print(paste("The estimated correlation length for mass in Land is:", pars_mass$thetam[2]*6371, "km", sep = " "))
print(paste("The estimated marginal variance for mass in Land is:", sqrt(pars_mass$thetam[1]), "mm/yr", sep = " "))

print(paste("The estimated correlation length for steric is:", pars_steric$rho_mode*6371, "km", sep = "  "))
print(paste("The estimated marginal variance for steric is:", sqrt(pars_steric$sigma_mode), "mm/yr", sep = "  "))



## ----predict, include=TRUE, echo = FALSE---------------------------------
steric_pred <- ress_2c_c$pred$steric
steric_pred$mean2 <- ifelse(abs(steric_pred$mean) > 19, sign(steric_pred$mean)*20, steric_pred$mean)
steric_pred$u2 <- ifelse(steric_pred$u > 3.8, 4, steric_pred$u)

ssh_pred <- ress_2c_c$pred$ssh
ssh_pred$mean2 <- ifelse(abs(ssh_pred$mean) > 19, sign(ssh_pred$mean)*20, ssh_pred$mean)
ssh_pred$u2 <- ifelse(ssh_pred$u > 3.8, 4, ssh_pred$u)

mass_pred <- ress_2c_c$pred$mass
mass_pred$mean2 <- ifelse(abs(mass_pred$mean) > 19, sign(mass_pred$mean)*20, mass_pred$mean)
mass_pred$u2 <- ifelse(mass_pred$u > 3.8, 4, mass_pred$u)

## plot the mass mean 
lattice::levelplot(mean2 ~ lon + lat, data = mass_pred, aspect = "iso", at = seq(-20, 20, 2),
                     panel = function(x,y,z,...){
                       lattice::panel.fill(col = "grey")
                       lattice::panel.levelplot(x,y,z,...)
                       map2 <- map("world2", interior = FALSE, plot = FALSE)
                       lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                     },
                     main = "The predicted mass trend (mm/yr ewh)", xlab = "longitude", ylab = "latitude")

## plot the mass uncertainty
lattice::levelplot(u2 ~ lon + lat, data = mass_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
                     panel = function(x,y,z,...){
                       lattice::panel.fill(col = "grey")
                       lattice::panel.levelplot(x,y,z,...)
                       map2 <- map("world2", interior = FALSE, plot = FALSE)
                       lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                     },
                     main = "The predited mass uncertainties (mm/yr ewh)", xlab = "longitude", ylab = "latitude")

## Plot the steric mean
lattice::levelplot(mean2 ~ lon + lat, data = steric_pred, aspect = "iso", at = seq(-20, 20, 2),
                     panel = function(x,y,z,...){
                       lattice::panel.fill(col = "grey")
                       lattice::panel.levelplot(x,y,z,...)
                       map2 <- map("world2", interior = FALSE, plot = FALSE)
                       lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                     },
                     main = "The predicted steric trend (mm/yr)", xlab = "longitude", ylab = "latitude")

## plot the steric uncertainty
lattice::levelplot(u2 ~ lon + lat, data = steric_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
                     panel = function(x,y,z,...){
                       lattice::panel.fill(col = "grey")
                       lattice::panel.levelplot(x,y,z,...)
                       map2 <- map("world2", interior = FALSE, plot = FALSE)
                       lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                     },
                     main = "The predited steric uncertainties (mm/yr)", xlab = "longitude", ylab = "latitude")

## plot the ssh mean 
lattice::levelplot(mean2 ~ lon + lat, data = ssh_pred, aspect = "iso", at = seq(-20, 20, 2),
                     panel = function(x,y,z,...){
                       lattice::panel.fill(col = "grey")
                       lattice::panel.levelplot(x,y,z,...)
                       map2 <- map("world2", interior = FALSE, plot = FALSE)
                       lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                     },
                     main = "The predicted SSH trend (mm/yr)", xlab = "longitude", ylab = "latitude")

## plot the ssh uncertainty
lattice::levelplot(u2 ~ lon + lat, data = ssh_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
                     panel = function(x,y,z,...){
                       lattice::panel.fill(col = "grey")
                       lattice::panel.levelplot(x,y,z,...)
                       map2 <- map("world2", interior = FALSE, plot = FALSE)
                       lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                     },
                     main = "The predited SSH uncertainties (mm/yr)", xlab = "longitude", ylab = "latitude")

## compare with grace
grace_df$predm <- ress_2c_c$pred$grace_pred$mean
grace_df$predu <- ress_2c_c$pred$grace_pred$u
grace_df$diff <- ress_2c_c$pred$grace_pred[,1] - grace_df$trendgia
grace_df$diff2  <- ifelse(abs(grace_df$diff) > 2, sign(grace_df$diff)*2.5, grace_df$diff )
grace_df_sp <- grace_df
coordinates(grace_df_sp) <- c("lon", "lat")
spplot(grace_df_sp, "diff2")
spplot(grace_df_sp, c("predu", "std"))
spplot(grace_df_sp, c("predm", "trendgia"))

grace_out <- grace_df[,c("lon", "lat", "predm", "predu", "trendgia")]
write.table(grace_out, file = "exp2_GRACE_output_c.txt", row.names = FALSE)

## ----local_constr--------------------------------------------------------
## Create the A matrix for the constraint
## Find which vertices fall in to South America
continents <- readOGR(dsn = paste0(dd,"WP1-BHM/maps/Continents"), layer = "continent")
SA <- continents[continents$CONTINENT == "South America",]
SAareas <- sapply(SA@polygons[[1]]@Polygons, function(x) x@area)
SA <- SA@polygons[[1]]@Polygons[SAareas > 10]
SA <- SpatialPolygons(list(Polygons(SA, ID = 1)))
proj4string(SA) <- proj4string(grace_df_sp) <- CRS(proj4string(grace_sp))
SA <- recenter(SA)

## ----globalmean, echo=FALSE----------------------------------------------
SA_id <- over(SA, grace_df_sp, returnList = TRUE)
eartharea <- sum(mass_pred$area)
oceanarea <-sum(subset(mass_pred, region == "Ocean")$area)
landarea <- sum(subset(mass_pred, region == "Land")$area)

predsum <- sum(mass_pred$mean*mass_pred$area)
predocean <- sum(subset(mass_pred, region == "Ocean")$mean*subset(mass_pred, region == "Ocean")$area)
predland <-sum(subset(mass_pred, region == "Land")$mean*subset(mass_pred, region == "Land")$area)
diffSA <- mean(SA_id[[1]]$diff)

## ----oceanmean, echo = FALSE---------------------------------------------
cat("GRACE predicted trend mean over the corresponding regions: \n",
    "\n",
    "Regions       |", "SA sum to zero",    "\n",
    "--------------------------------------- \n",
    "Global        | ", predsum/eartharea,    "\n",
    "Ocean         | ", predocean/oceanarea,     "\n",
    "Land          |", predland/landarea,        "\n",
    "South America |", diffSA,          "\n")


## ----2ncdf, include = FALSE, eval = FALSE--------------------------------
## source(paste0(wd,"gmrcode/BHM_sphere/ncdfio.R"))
## mass_pred <- mass_pred[order(mass_pred$lat, mass_pred$lon),]
## 
## ## save mass
## df2ncdf(df = mass_pred[, c("lon", "lat", "mean", "u")], fname = "exp2c_c_mass.nc",
##         vars = list(c("mean", "trend", "BHM predicted mass trend", "mm/yr"),
##              c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
##         title = "BHM-predicted mass",
##         append = FALSE)
## 
## ## save steric
## df2ncdf(df = ress_2c_c$pred$steric, fname = "exp2c_c_steric.nc",
##         vars = list(c("mean", "trend", "BHM predicted steric trend", "mm/yr"),
##              c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
##         title = "BHM-predicted steric",
##         append = FALSE)
## 
## ## save SSH
## df2ncdf(df = ress_2c_c$pred$ssh, fname = "exp2c_c_ssh.nc",
##         vars = list(c("mean", "trend", "BHM predicted SSH trend", "mm/yr"),
##              c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
##         title = "BHM-predicted SSH",
##         append = FALSE)

