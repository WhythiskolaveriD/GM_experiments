## ----load_pack, message = FALSE, warning=FALSE---------------------------
library(sp); library(INLA); library(GEOmap)
library(ggplot2); library(grid); library(gridExtra)
source("functions.R")

## ----load_data0, include=FALSE, eval = TRUE, cache = TRUE----------------
#### 1 Load GIA prior
if(Sys.info()["sysname"] == "Windows"){
  ice6g <- read.table("Z:/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)
}else if(Sys.info()["sysname"] == "Ubuntu"){
  ice6g <- read.table("~/globalmass/WP2-SolidEarth/GIAforwardModels/textfiles/GIA_Pel-6-VM5.txt", header = T)
}else{
  ice6g <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)
}

polycoords <- ice6g[,c(6:13, 6,7)] 
plist <- lapply(ice6g$ID, 
                function(x) Polygons(list(Polygon(cbind(lon = as.numeric(polycoords[x, c(1,3,5,7,9)]), 
                                                        lat = as.numeric(polycoords[x, c(2,4,6,8,10)])))), ID = x))
Plist <- SpatialPolygons(plist, proj4string = CRS("+proj=longlat"))


#### 2 Load GPS data
if(Sys.info()["sysname"] == "Windows"){
  GPSV4b <- read.table("Z:/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)
}else if(Sys.info()["sysname"] == "Ubuntu"){
  GPSV4b <- read.table("~/globalmass/WP2-SolidEarth/GPS/NGL/BHMinputFiles/GPS_v04b.txt", header = T)
}else{
  GPSV4b <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)
}

## ----gia_data, include=TRUE, eval = FALSE--------------------------------
## ice6g <- read.table("ice6g.txt", header = T)
## polycoords <- ice6g[,c(6:13, 6,7)]
## plist <- lapply(ice6g$ID,
##                 function(x) Polygons(list(Polygon(cbind(lon = as.numeric(polycoords[x, c(1,3,5,7,9)]),
##                                                         lat = as.numeric(polycoords[x, c(2,4,6,8,10)])))), ID = x))
## Plist <- SpatialPolygons(plist, proj4string = CRS("+proj=longlat"))
## ## Note that in this data set the grid boundaries are defined on [-0.5, - 359.5] and the centre points on [0, 359].
## ## Corresponding transformation is needed in the following in geometery operations.

## ----fibo_mesh, include=TRUE---------------------------------------------
fibo_points <- fiboSphere(N = 3e4, LL = FALSE)
mesh <- inla.mesh.2d(loc = fibo_points, cutoff = 0.01, max.edge = 0.5)
summary(mesh)

## ----GIA_mesh, include=TRUE, eval = TRUE, cache = TRUE-------------------
meshLL <- Lxyz2ll(list(x=mesh$loc[,1], y = mesh$loc[,2], z = mesh$loc[,3]))
meshLL$lon <- ifelse(meshLL$lon >= -0.5, meshLL$lon,meshLL$lon + 360)
mesh_sp <- SpatialPoints(data.frame(lon = meshLL$lon, lat = meshLL$lat), proj4string = CRS("+proj=longlat")) 
mesh_idx <- over(mesh_sp, Plist)
GIA_prior <- ice6g$trend[mesh_idx]

## ----loadGPS, include = TRUE, eval = FALSE-------------------------------
## GPSV4b <- read.table("GPS_Versionx.txt", header = T)
## ## The original GPS data use lon in [0, 360]

## ----GPSdetrend, include = TRUE------------------------------------------
GPS_data <- GPSV4b
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_data$lat, lon = GPS_data$lon))
GPS_sp <- SpatialPoints(data.frame(lon = ifelse(GPS_data$lon>359.5, GPS_data$lon - 360, GPS_data$lon), 
                                   lat = GPS_data$lat), proj4string = CRS("+proj=longlat"))
GPS_idx <- over(GPS_sp, Plist)
GPS_mu <- ice6g$trend[GPS_idx]
GPS_data$trend0 <- GPS_data$trend - GPS_mu

## ----prior---------------------------------------------------------------
## Priors mean and variance for the parameters: rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

## Transform the parameters for the SPDE_GMRF approximation
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

## ----inla_setup, include = TRUE, cache = TRUE----------------------------
Mesh_GIA <- mesh
## The SPDE model
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0
GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

## Link the process to observations and predictions
A_data <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
A_pred <- inla.spde.make.A(mesh = Mesh_GIA, loc = rbind(GPS_loc, Mesh_GIA$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=GPS_data$trend0), A = list(A_data),
                     effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

## Fix the GPS errors
hyper <- list(prec = list(fixed = TRUE, initial = 0))
formula = y ~ -1 +  f(GIA, model = GIA_spde)
prec_scale <- c(1/GPS_data$std^2, rep(1, nrow(A_pred)))

##----inla_run, include = TRUE, eval = FALSE------------------------------
## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))
save(res_inla, file = "/./projects/GlobalMass/WP1-BHM/Experiment1b/GIA_RGL/res1.RData")

INLA_pred <- res_inla$summary.linear.predictor

## ----inla_res, include = TRUE, cache=TRUE--------------------------------
INLA_pred <- res_inla$summary.linear.predictor

## Extract and project predictions
pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
GPS_idx <- pred_idx[1:nrow(GPS_data)]
GIA_idx <- pred_idx[-(1:nrow(GPS_data))]

## GPS 
GPS_u <- INLA_pred$sd[GPS_idx]
GPS_pred <- data.frame(lon = GPS_data$lon, lat = GPS_data$lat, u = GPS_u)

## GIA
GIA_diff <- INLA_pred$mean[GIA_idx] 
GIA_m <- GIA_diff + GIA_prior
GIA_u <- INLA_pred$sd[GIA_idx]
proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(360,180), xlim = c(0,360), ylim = c(-90, 90))
GIA_grid <- expand.grid(proj$x, proj$y)
GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                       diff = as.vector(inla.mesh.project(proj, as.vector(GIA_diff))),
                       mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                       u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))))

# inla predicts all zero
ress <- list(res_inla = res_inla, spde = GIA_spde, st = stGIA, 
            mesh = Mesh_GIA, GPS_pred = GPS_pred, GIA_pred = GIA_pred)

## ----hyper, include=TRUE-------------------------------------------------
res_inla <- ress$res_inla
GIA_spde <- ress$spde
pars_GIA <- inla.spde2.result(res_inla, "GIA", GIA_spde, do.transf=TRUE)
theta_mean <- pars_GIA$summary.theta$mean
theta_sd <- pars_GIA$summary.theta$sd

## Find the mode of rho and sigma^2
lrho_mode <- pars_GIA$summary.log.range.nominal$mode
lrho_mean <- pars_GIA$summary.log.range.nominal$mean
lrho_sd <- pars_GIA$summary.log.range.nominal$sd
rho_mode <- exp(lrho_mean - lrho_sd^2)

lsigma_mode <- pars_GIA$summary.log.variance.nominal$mode
lsigma_mean <- pars_GIA$summary.log.variance.nominal$mean
lsigma_sd <- pars_GIA$summary.log.variance.nominal$sd
sigma_mode <- exp(lsigma_mean - lsigma_sd^2)

plot(pars_GIA$marginals.range.nominal[[1]], type = "l",
     main = bquote(bold(rho("mode") == .(round(rho_mode, 4))))) # The posterior from inla output
plot(pars_GIA$marginals.variance.nominal[[1]], type = "l", 
     main = bquote(bold({sigma^2}("mode") == .(round(sigma_mode, 4))))) # The posterior from inla output

## The estimated correlation length is about 568km
rho_mode*6371

## ----predict, include=TRUE-----------------------------------------------
GPS_pred <- ress$GPS_pred
GIA_pred <- ress$GIA_pred

map_prior <- map_res(data = ice6g, xname = "x_center", yname = "y_center", fillvar = "trend", 
                    limits = c(-7, 22), title = "Prior GIA mean field")

## Plot the GIA predicted mean
map_GIA <- map_res(data = GIA_pred, xname = "lon", yname = "lat", fillvar = "mean", 
                  limits = c(-7, 22), title = "Predicted GIA")

## Plot the GIA difference map
map_diff <- map_res(data = GIA_pred, xname = "lon", yname = "lat", fillvar = "diff", 
                  limits = c(-8, 8), title = "GIA difference: Updated - Prior")

## Plot the GIA difference map
map_sd <- map_res(data = GIA_pred, xname = "lon", yname = "lat", fillvar = "u", 
                    colpal = colorRamps::matlab.like(12), title = "Predicted uncertainties")

 
## Display
print(map_prior)
print(map_GIA)
print(map_diff)
print(map_sd)

