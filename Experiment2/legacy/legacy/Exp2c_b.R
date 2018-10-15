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
  ## any server  
  wd <- "~/"
  dd <- "/./projects/GlobalMass/"
}



library(rgdal); library(sp);library(GEOmap)
library(INLA)
library(ncdf4)
source(paste0(wd,"gmrcode/BHM_sphere/functions.R"))
source(paste0(wd, "gmrcode/BHM_sphere/partition_fun.R"))
source(paste0(wd, "gmrcode/BHM_sphere/pp_functions.R"))

## ----loaddata, include = FALSE-------------------------------------------
## The altimetry data
alt_nc <- nc_open(paste0(dd,"WP3-Ocean/BHMinputs/SSH/trend_SSH_CCI_200501_201512.nc"))
lon <- ncvar_get(alt_nc, "lon")
lat <- ncvar_get(alt_nc, "lat")
trend_ssh <- ncvar_get(alt_nc, "trend_ssh") #note that there are NAs for land datat
err_ssh <- ncvar_get(alt_nc, "err")
alt_data <- data.frame(trend_ssh = as.numeric(trend_ssh), err_ssh = as.numeric(err_ssh))
alt_data$lon <- rep(lon, 180) 
alt_data$lat <- rep(lat, each = 360)
alt_data2 <- na.omit(alt_data)
alt_loc <- do.call(cbind, Lll2xyz(lon = alt_data2$lon, lat = alt_data2$lat))

## Load BDV's GRACE data -- not the data are not evenly spaced!
grace_nc <- nc_open(paste0(dd, "WP2-SolidEarth/Bramha/BHMinput/GRACE_trends_Oceanadded_Bramha.nc"))
print(grace_nc)
lon <- ncvar_get(grace_nc, "Long")
lat <- ncvar_get(grace_nc, "Lat")
trend_grace <- ncvar_get(grace_nc, "trend") #note that there are NAs for land datat
err_grace <- ncvar_get(grace_nc, "uncert")

grace_df <- data.frame(trend = as.numeric(trend_grace), std = as.numeric(err_grace), lon= lon, lat = lat)
grace_bdv <- grace_df
coordinates(grace_bdv) <- c("lon", "lat")

## Load BDV's ICE6G-VM5 in ewh mm/yr
gia_ewh <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA_BDV/GIA_ewh_ICE6G"), sep = ",", head = TRUE)
gia_ewh$lat <- -gia_ewh$lat
names(gia_ewh) <- c("lat", "lon", "trend")
gia_ewh$trend <- gia_ewh$trend*1000


## Load the Land polygon (-180, 180)
Land <- readOGR(dsn = paste0(dd, "WP1-BHM/maps/ne_110m_land"), layer = "ne_110m_land")
## Remove small islands in the low lattidure(lat < 70) regrions (area < 300*300km)
areas <- sapply(Land@polygons, slot, "area")
lats <- sapply(Land@polygons, function(x){slot(x, "labpt")[2]})
smallareas <- areas <= 10
lowlats <- lats <= 70
lowsmall <- lowlats & smallareas
newids <- which(!lowsmall)
newLand <- SpatialPolygons(Land@polygons[newids], proj4string = CRS(proj4string(Land)))

## The ICE6G-VM5 in vlm mm/yr
ice6g_vlm <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt"), header = T)
names(ice6g_vlm)[2:3] <- c("lon", "lat")
ice6g_vlm$lon <- ifelse(ice6g_vlm$lon >180, ice6g_vlm$lon - 360, ice6g_vlm$lon)
coordinates(ice6g_vlm) <-c("lon", "lat")
gridded(ice6g_vlm) <- TRUE
proj4string(ice6g_vlm) <- CRS(proj4string(Land))
ice6g_vlm$areas <- geosphere::areaPolygon(as(ice6g_vlm, "SpatialPolygons"))/(1000^2) # calculated the grid areas
Lid <- unlist(over(newLand, as(ice6g_vlm, "SpatialPoints"), returnList = T)) # Find the grid in land
vlm <- sum(ice6g_vlm$trend[-Lid]*ice6g_vlm$areas[-Lid]/sum(ice6g_vlm$areas[-Lid])) # the average vlm in the ocean
print(paste("The vlm adjustment is ", vlm, "mm/yr"))

## Adjust the altimetry data by vlm
alt_data2$trend_vlm <- alt_data2$trend_ssh - vlm

## ----grace_gia-----------------------------------------------------------
## GRACE sum to zero -- equal area -- simple average
mean(grace_bdv$trend)

## GIA_ewh sum to zero
gia_sp <- gia_ewh
coordinates(gia_sp) <- c("lon", "lat")
gridded(gia_sp) <- TRUE
gia_sp$areas <- geosphere::areaPolygon(as(gia_sp, "SpatialPolygons"))/(1000^2)
sum(gia_sp$trend*gia_sp$areas)/sum(gia_sp$areas)

## Remove GIA from GRACE -- use grace - gia grid value which grace falls in
proj4string(gia_sp)<- proj4string(grace_bdv) <- CRS("+proj=longlat")
gridded(gia_sp) <- TRUE
giaG <- over(grace_bdv, as(gia_sp, "SpatialPolygons"))
grace_df$trendgia <- grace_bdv$trendgia <- grace_bdv$trend - gia_sp$trend[giaG]
mean(grace_bdv$trendgia)

## ----learnmass-----------------------------------------------------------
grace_coords <- grace_bdv@coords
grace_coords[,1] <- ifelse(grace_coords[,1]> 180, grace_coords[,1]-360, grace_coords[,1])
glid <- unlist(over(newLand, SpatialPoints(coords= grace_coords, proj4string = CRS(proj4string(newLand))), returnList = TRUE))
goid <- (1:nrow(grace_coords))[-glid]
# v_ocean <- gstat::variogram(trendgia ~ 1, data = grace_bdv[goid,])
# v_land <- gstat::variogram(trendgia ~ 1, data = grace_bdv[-goid,])
# plot(v_ocean)
# plot(v_land)
# -- the range for ocean is likely to be about 1200km 
# -- the variance is likely to be around (9mm/yr)^2

# -- the range for land is likely to be 600km
# -- the variance is (19mm/yr)^2


## ----prior---------------------------------------------------------------
## Priors for mass
# -- for the ocean
mu_rg1 <- 1200/6371 # mean of rho
v_rg1 <- 0.1 # vague variace for rho
mu_sg1 <- 9 # mean of sigma
v_sg1 <- 18^2 # vaue variance for sigma
# -- for the land
mu_rg2 <- 300/6371 # mean of rho
v_rg2 <- 0.1 # vague variace for rho
mu_sg2 <- 19 # mean of sigma
v_sg2 <- 38^2 # vaue variance for sigma

## Priors for steric
mu_r2 <- 2000/6371 # mean of rho
v_r2 <- 1 # vague variace for rho
mu_s2 <- 15 # mean of sigma
v_s2 <- 30^2 # vaue variance for sigma

## Transform the parameters for a lognormal distribution
## -- mass
trhog1 <- Tlognorm(mu_rg1, v_rg1)
tsigmag1 <- Tlognorm(mu_sg1, v_sg1)
trhog2 <- Tlognorm(mu_rg2, v_rg2)
tsigmag2 <- Tlognorm(mu_sg2, v_sg2)

lsigmag1 <- tsigmag1[1]
lrhog1 <- trhog1[1]


lsigmag2 <- tsigmag2[1]
lrhog2 <- trhog2[1]


## -- steric
trho2 <- Tlognorm(mu_r2, v_r2)
tsigma2 <- Tlognorm(mu_s2, v_s2)

lsigma2 <- tsigma2[1]
theta2_ss <- tsigma2[2]
lrho2 <- trho2[1]
theta2_rs <- trho2[2]
lkappa2 <- log(8)/2 - lrho2
ltau2 <- 0.5*log(1/(4*pi)) - lsigma2 - lkappa2

## ----ssh_mesh------------------------------------------------------------
## Genereate Fibonacci points on the sphere
fibo_points <- fiboSphere(N = 3e4, LL = FALSE)
mesh0 <- inla.mesh.2d(loc = fibo_points, cutoff = 0.01, max.edge = 1)
summary(mesh0)

## ----ssh_mesh2-----------------------------------------------------------
mesh0 <- dt.mesh.addon.posTri(mesh = mesh0, globe = TRUE) 
Tlonlat <- Lxyz2ll(list(x = mesh0$posTri[,1], y = mesh0$posTri[,2], z = mesh0$posTri[,3]))
mesh0$Trill <- cbind(lon = Tlonlat$lon, lat =Tlonlat$lat)
TinLand <- unlist(over(newLand, SpatialPoints(coords=mesh0$Trill, proj4string = CRS(proj4string(newLand))), returnList=T))
TAll <- 1:mesh0$t
TinOcean <- TAll[-TinLand]
Omega = dt.Omega(list(TinOcean, 1:mesh0$t), mesh0)

mesh_ocean <- mesh.sub(mesh0, Omega, 1)
mesh_land <- mesh.sub(mesh0, Omega, 2)
summary(mesh_ocean)
summary(mesh_land)

## ----spde----------------------------------------------------------------
## mass
massQ <- pp.create.Q(meshes = list(mesh_ocean, mesh_land), initial.theta = c(lrhog1, lrhog2, lsigmag1, lsigmag2))
## Transform the parameters for the SPDE_GMRF approximation
prior <- list(sigma = rbind(tsigmag1, tsigmag2), rho = rbind(trhog1, trhog2))
log.prior <- pp.create.prior.log.norm(prior.param = prior) 
mass_spde <- pp.inla.model(Q = massQ, log.prior=log.prior)

## steric
steric_spde <- inla.spde2.matern(mesh_ocean, B.tau = matrix(c(ltau2, -1, 1),1,3), B.kappa = matrix(c(lkappa2, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta2_ss), sqrt(1/theta2_rs)))

## ----link_data-----------------------------------------------------------
## The data
ygrace_Ocean <- grace_bdv$trendgia[goid]
ygrace_Land <- grace_bdv$trendgia[-goid]
yaltv<- alt_data2$trend_vlm

## The A matrice for linking observation and processes

## Link GRACE to mass
grace_locO <- do.call(cbind, Lll2xyz(lat = grace_bdv@coords[goid,2], lon = grace_bdv@coords[goid,1]))
grace_locL <- do.call(cbind, Lll2xyz(lat = grace_bdv@coords[-goid,2], lon = grace_bdv@coords[-goid,1]))
A_graceO <- inla.spde.make.A(mesh = mesh_ocean, loc = grace_locO)
A_graceL <- inla.spde.make.A(mesh = mesh_land, loc = grace_locL)
A_grace_m <- bdiag(A_graceO, A_graceL)

## Link Altimetry and grace to steric and mass
A_altv_s <- inla.spde.make.A(mesh = mesh_ocean, loc = alt_loc)
A_altv_m <- inla.spde.make.A(mesh = mesh_ocean, loc = alt_loc)
A_altv_m <- cBind(A_altv_m, Matrix(0, nrow = nrow(A_altv_m), ncol = ncol(A_grace_m) - ncol(A_altv_m)))

## The A matrices for predict ssh
A_ssh_s <- Diagonal(n= steric_spde$n.spde, x =1)
A_ssh_m <- cBind(Diagonal(n= mesh_ocean$n, x =1), Matrix(0, nrow = mesh_ocean$n, ncol = mesh_land$n))

## The errors
prec_scale <- c(1/grace_bdv$std[goid]^2, 1/grace_bdv$std[-goid]^2, 1/alt_data2$err_ssh^2, rep(1, nrow(A_ssh_s)))

## Create the mass and steric stacks
stkmass <- inla.stack(data = list(y=c(ygrace_Ocean, ygrace_Land)), A = list(A_grace_m),
                     effects = list(mass = 1:(mesh_ocean$n + mesh_land$n)), 
                      remove.unused = FALSE, tag = "mass")

stksteric <- inla.stack(data = list(y=yaltv), A = list(A_altv_m, A_altv_s),
                     effects = list(list(mass = 1:(mesh_ocean$n + mesh_land$n)),
                                    list(steric = 1:steric_spde$n.spde)), 
                     remove.unused = FALSE, tag = "steric")

## Create the stack for predicting SSH
stkssh <- inla.stack(data = list(y=NA), A = list(A_ssh_m, A_ssh_s),
                     effects = list(list(mass = 1:(mesh_ocean$n + mesh_land$n)), 
                                    list(steric = 1:steric_spde$n.spde)), 
                     remove.unused = FALSE, tag = "ssh")
## join the two stacks
stkall <- inla.stack(stkmass,  stksteric, stkssh, remove.unused = FALSE)

## ----inla_run, include = TRUE, eval = FALSE------------------------------
## ## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))

## The formular -- we add the constraint that mass change sum to zero
## constraint 1 -- vertices sum to zero
A1 <- matrix(1, nrow = 1, ncol = (mesh_ocean$n + mesh_land$n))
## Constraint 2 -- grace loc sum to zero
A2 <- matrix(colSums(A_grace_m), nrow = 1)
## Constraint 3 -- longlat grid sum to zero
gridll <- ice6g_vlm@coords
gridxyz <- do.call(cbind, Lll2xyz(lat = gridll[,2], lon = gridll[,1]))
gridxyzO <- gridxyz[-Lid,]
gridxyzL <- gridxyz[Lid,]
A_gridO <- inla.spde.make.A(mesh = mesh_ocean, loc = gridxyzO)
A_gridL <- inla.spde.make.A(mesh = mesh_land, loc = gridxyzL)
A_grid_m <- bdiag(A_gridO, A_gridL)
weights <- ice6g_vlm$areas/sum(ice6g_vlm$areas)
A3 <- matrix(weights, nrow = 1) %*% A_grid_m
A <- as.matrix(rbind(A1, A2, A3))
## The formular
formula = y ~ -1 + f(mass, model = mass_spde,
                     extraconstr = list(A = A, e = c(0, 0, 0))) + f(steric, model = steric_spde)

## Run INLA

res_inla <- inla(formula, data = inla.stack.data(stkall), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stkall), compute =TRUE), verbose = TRUE)

## ----inla_res, include = TRUE, eval = FALSE, echo = FALSE----------------
## idices of the predicted processes
mass_idx <- inla.stack.index(stkall, tag = "mass")$effects
mass_idxg <- inla.stack.index(stkall, tag = "mass")$data
steric_idx <- inla.stack.index(stkall, tag = "steric")$effects[-c(1:(mesh_ocean$n + mesh_land$n))]
ssh_idx <- inla.stack.index(stkall, tag = "ssh")$data


## Extract the predictions
INLA_mass <- res_inla$summary.random$mass
INLA_steric <- res_inla$summary.random$steric
INLA_ssh <- res_inla$summary.linear.predictor[ssh_idx,]
INLA_grace <- res_inla$summary.linear.predictor[mass_idxg,]

## mass
## mass Ocean
massO_mean <- INLA_mass$mean[1:mesh_ocean$n]
massO_u <- INLA_mass$sd[1:mesh_ocean$n]
proj <- inla.mesh.projector(mesh_ocean, projection = "longlat", dims = c(360,180), xlim = c(0.5,359.5), ylim = c(-89.5, 89.5))
mass_grid <- expand.grid(proj$x, proj$y)
massO_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                       mean = as.vector(inla.mesh.project(proj, as.vector(massO_mean))),
                       u = as.vector(inla.mesh.project(proj, as.vector(massO_u))),
                       region = "Ocean")

coords <- cbind(massO_pred$lon, massO_pred$lat)
mass_sp <- SpatialPoints(coords=coords, proj4string = CRS("+proj=longlat"))
mass_spg <- points2grid(mass_sp)
grid_areas <- geosphere::areaPolygon(as(mass_spg, "SpatialPolygons"))/(1000^2)
massO_pred$area <- grid_areas

## mass Land
massL_mean <- INLA_mass$mean[-(1:mesh_ocean$n)]
massL_u <- INLA_mass$sd[-(1:mesh_ocean$n)]
proj <- inla.mesh.projector(mesh_land, projection = "longlat", dims = c(360,180), xlim = c(0.5,359.5), ylim = c(-89.5, 89.5))
mass_grid <- expand.grid(proj$x, proj$y)
massL_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                       mean = as.vector(inla.mesh.project(proj, as.vector(massL_mean))),
                       u = as.vector(inla.mesh.project(proj, as.vector(massL_u))),
                       region = "Land")
massL_pred$area <- grid_areas
## mass all
mass_pred <- na.omit(rbind(massO_pred, massL_pred))

## mass at GRACE resolution
grace_pred <- data.frame(mean = INLA_grace$mean, u = INLA_grace$sd)

## steric
steric_mean <- INLA_steric$mean
steric_u <- INLA_steric$sd
proj2 <- inla.mesh.projector(mesh_ocean, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-89.5, 89.5))
steric_grid <- expand.grid(proj2$x, proj2$y)
steric_pred <- data.frame(lon = steric_grid[,1], lat = steric_grid[,2],
                       mean = as.vector(inla.mesh.project(proj2, as.vector(steric_mean))),
                       u = as.vector(inla.mesh.project(proj2, as.vector(steric_u))))

## ssh
ssh_mean <- INLA_ssh$mean + vlm
ssh_u <- INLA_ssh$sd
ssh_pred <- data.frame(lon = steric_grid[,1], lat = steric_grid[,2],
                       mean = as.vector(inla.mesh.project(proj2, as.vector(ssh_mean))),
                       u = as.vector(inla.mesh.project(proj2, as.vector(ssh_u))))

ress_2c_b <- list(res_inla = res_inla, st = stkall,
                spde = list(mass = mass_spde, steric = steric_spde),
                mesh = list(mass = mesh0, steric = mesh_ocean),
                pred = list(mass = mass_pred, steric = steric_pred, ssh = ssh_pred, grace_pred = grace_pred))

save(ress_2c_b, file = "/./projects/GlobalMass/WP1-BHM/Experiment2c/Exp2c_b.RData")


## ----hyper, include=TRUE, echo = FALSE-----------------------------------
# load(paste0(dd, "WP1-BHM/Experiment2c/Exp2c_b.RData"))
# 
# pars_mass <- marginal_par(res = ress_2c_b, process = "mass", SPDE2 = FALSE, 
#                          theta.names = c("rho_O","rho_L", "sigma_O", "sigma_L"), plot = TRUE)
# pars_steric <- marginal_par(res = ress_2c_b, process = "steric", SPDE2 = TRUE, 
#                             theta.names = c("rho_S","sigma_S"), plot = TRUE)
# 
# ## The posterior modes
# print(paste("The estimated correlation length for mass in Ocean is:", pars_mass$thetam[1]*6371, "km", sep = " "))
# print(paste("The estimated marginal variance for mass in Ocean is:", sqrt(pars_mass$thetam[3]), "mm/yr", sep = " "))
# 
# print(paste("The estimated correlation length for mass in Land is:", pars_mass$thetam[2]*6371, "km", sep = " "))
# print(paste("The estimated marginal variance for mass in Land is:", sqrt(pars_mass$thetam[4]), "mm/yr", sep = " "))
# 
# print(paste("The estimated correlation length for steric is:", pars_steric$rho_mode*6371, "km", sep = "  "))
# print(paste("The estimated marginal variance for steric is:", sqrt(pars_steric$sigma_mode), "mm/yr", sep = "  "))



## ----predict, include=TRUE, echo = FALSE---------------------------------
steric_pred <- ress_2c_b$pred$steric
steric_pred$mean2 <- ifelse(abs(steric_pred$mean) > 19, sign(steric_pred$mean)*20, steric_pred$mean)
steric_pred$u2 <- ifelse(steric_pred$u > 3.8, 4, steric_pred$u)

ssh_pred <- ress_2c_b$pred$ssh
ssh_pred$mean2 <- ifelse(abs(ssh_pred$mean) > 19, sign(ssh_pred$mean)*20, ssh_pred$mean)
ssh_pred$u2 <- ifelse(ssh_pred$u > 3.8, 4, ssh_pred$u)

mass_pred <- ress_2c_b$pred$mass
mass_pred$mean2 <- ifelse(abs(mass_pred$mean) > 19, sign(mass_pred$mean)*20, mass_pred$mean)
mass_pred$u2 <- ifelse(mass_pred$u > 3.8, 4, mass_pred$u)

# ## plot the mass mean 
# lattice::levelplot(mean2 ~ lon + lat, data = mass_pred, aspect = "iso", at = seq(-20, 20, 2),
#                      panel = function(x,y,z,...){
#                        lattice::panel.fill(col = "grey")
#                        lattice::panel.levelplot(x,y,z,...)
#                        map2 <- map("world2", interior = FALSE, plot = FALSE)
#                        lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
#                      },
#                      main = "The predicted mass trend (mm/yr ewh)", xlab = "longitude", ylab = "latitude")
# 
# ## plot the mass uncertainty
# lattice::levelplot(u2 ~ lon + lat, data = mass_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
#                      panel = function(x,y,z,...){
#                        lattice::panel.fill(col = "grey")
#                        lattice::panel.levelplot(x,y,z,...)
#                        map2 <- map("world2", interior = FALSE, plot = FALSE)
#                        lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
#                      },
#                      main = "The predited mass uncertainties (mm/yr ewh)", xlab = "longitude", ylab = "latitude")
# 
# ## Plot the steric mean
# lattice::levelplot(mean2 ~ lon + lat, data = steric_pred, aspect = "iso", at = seq(-20, 20, 2),
#                      panel = function(x,y,z,...){
#                        lattice::panel.fill(col = "grey")
#                        lattice::panel.levelplot(x,y,z,...)
#                        map2 <- map("world2", interior = FALSE, plot = FALSE)
#                        lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
#                      },
#                      main = "The predicted steric trend (mm/yr)", xlab = "longitude", ylab = "latitude")
# 
# ## plot the steric uncertainty
# lattice::levelplot(u2 ~ lon + lat, data = steric_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
#                      panel = function(x,y,z,...){
#                        lattice::panel.fill(col = "grey")
#                        lattice::panel.levelplot(x,y,z,...)
#                        map2 <- map("world2", interior = FALSE, plot = FALSE)
#                        lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
#                      },
#                      main = "The predited steric uncertainties (mm/yr)", xlab = "longitude", ylab = "latitude")
# 
# ## plot the ssh mean 
# lattice::levelplot(mean2 ~ lon + lat, data = ssh_pred, aspect = "iso", at = seq(-20, 20, 2),
#                      panel = function(x,y,z,...){
#                        lattice::panel.fill(col = "grey")
#                        lattice::panel.levelplot(x,y,z,...)
#                        map2 <- map("world2", interior = FALSE, plot = FALSE)
#                        lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
#                      },
#                      main = "The predicted SSH trend (mm/yr)", xlab = "longitude", ylab = "latitude")
# 
# ## plot the ssh uncertainty
# lattice::levelplot(u2 ~ lon + lat, data = ssh_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
#                      panel = function(x,y,z,...){
#                        lattice::panel.fill(col = "grey")
#                        lattice::panel.levelplot(x,y,z,...)
#                        map2 <- map("world2", interior = FALSE, plot = FALSE)
#                        lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
#                      },
#                      main = "The predited SSH uncertainties (mm/yr)", xlab = "longitude", ylab = "latitude")

## compare with grace
grace_df2 <- rbind(grace_df[goid,], grace_df[-goid,])
grace_df2$predm <- ress_2c_b$pred$grace_pred$mean
grace_df2$predu <- ress_2c_b$pred$grace_pred$u
grace_df2$diff <- ress_2c_b$pred$grace_pred[,1] - grace_df2$trendgia
grace_df2$diff2  <- ifelse(abs(grace_df2$diff) > 2, sign(grace_df2$diff)*2.5, grace_df2$diff )
# grace_df2_sp <- grace_df2
# coordinates(grace_df2_sp) <- c("lon", "lat")
# spplot(grace_df2_sp, "diff2")
# spplot(grace_df2_sp, c("predu", "std"))
# spplot(grace_df2_sp, c("predm", "trendgia"))

grace_out <- grace_df2[,c("lon", "lat", "predm", "predu", "trendgia")]
write.table(grace_out, file = "exp2_GRACE_output_b.txt", row.names = FALSE)

## ----local_constr--------------------------------------------------------
## Create the A matrix for the constraint
## Find which vertices fall in to South America
# continents <- readOGR(dsn = paste0(dd,"WP1-BHM/maps/Continents"), layer = "continent")
# SA <- continents[continents$CONTINENT == "South America",]
# SAareas <- sapply(SA@polygons[[1]]@Polygons, function(x) x@area)
# SA <- SA@polygons[[1]]@Polygons[SAareas > 10]
# SA <- SpatialPolygons(list(Polygons(SA, ID = 1)))
# proj4string(SA) <- proj4string(grace_df2_sp) <- CRS(proj4string(grace_bdv))
# SA <- recenter(SA)
# 
# ## ----globalmean, echo=FALSE----------------------------------------------
# SA_id <- over(SA, grace_df2_sp, returnList = TRUE)
# eartharea <- sum(mass_pred$area)
# oceanarea <-sum(subset(mass_pred, region == "Ocean")$area)
# landarea <- sum(subset(mass_pred, region == "Land")$area)
# 
# predsum <- sum(mass_pred$mean*mass_pred$area)
# predocean <- sum(subset(mass_pred, region == "Ocean")$mean*subset(mass_pred, region == "Ocean")$area)
# predland <-sum(subset(mass_pred, region == "Land")$mean*subset(mass_pred, region == "Land")$area)
# diffSA <- mean(SA_id[[1]]$diff)
# 
# ## ----oceanmean, echo = FALSE---------------------------------------------
# cat("GRACE predicted trend mean over the corresponding regions: \n",
#     "\n",
#     "Regions       |", "SA sum to zero",    "\n",
#     "--------------------------------------- \n",
#     "Global        | ", predsum/eartharea,    "\n",
#     "Ocean         | ", predocean/oceanarea,     "\n",
#     "Land          |", predland/landarea,        "\n",
#     "South America |", diffSA,          "\n")
# 

## ----2ncdf, include = FALSE, eval = FALSE--------------------------------
source(paste0(wd,"gmrcode/BHM_sphere/ncdfio.R"))
mass_pred <- mass_pred[order(mass_pred$lat, mass_pred$lon),]

## save mass
df2ncdf(df = mass_pred[, c("lon", "lat", "mean", "u")], fname = "exp2c_b_mass.nc",
        vars = list(c("mean", "trend", "BHM predicted mass trend", "mm/yr"),
             c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
        title = "BHM-predicted mass",
        append = FALSE)

## save steric
df2ncdf(df = ress_2c_b$pred$steric, fname = "exp2c_b_steric.nc",
        vars = list(c("mean", "trend", "BHM predicted steric trend", "mm/yr"),
             c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
        title = "BHM-predicted steric",
        append = FALSE)

## save SSH
df2ncdf(df = ress_2c_b$pred$ssh, fname = "exp2c_b_ssh.nc",
        vars = list(c("mean", "trend", "BHM predicted SSH trend", "mm/yr"),
             c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
        title = "BHM-predicted SSH",
        append = FALSE)

