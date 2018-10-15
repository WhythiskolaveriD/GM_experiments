## Load packages
library(sp); library(INLA); library(GEOmap)
source("glbm/BHM_sphere/functions.R")
source("glbm/BHM_sphere/partition_fun.R")

############################################################################
## 0 Data and Mesh
############################################################################

## 1 Load GIA prior
ice6g <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)
polycoords <- ice6g[,c(6:13, 6,7)] 
plist <- lapply(ice6g$ID, 
                function(x) Polygons(list(Polygon(cbind(lon = as.numeric(polycoords[x, c(1,3,5,7,9)]), 
                                                        lat = as.numeric(polycoords[x, c(2,4,6,8,10)])))), ID = x))
Plist <- SpatialPolygons(plist, proj4string = CRS("+proj=longlat"))
## Note that in this data set the grid boundaries are defined on [-0.5, - 359.5] and the centre points on [0, 359].
## Corresponding transformation is needed in the following in geometery operations.

## 2 Load GPS data
GPSV4b <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)


## 3 Priors mean and variance for the parameters: rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2
## Transform the parameters for the SPDE_GMRF approximation
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

## 4 Generate mesh
fibo_points <- fiboSphere(N = 12960, L0 = TRUE)
fibo_points_xyz <- do.call(cbind, Lll2xyz(lat = fibo_points[,2], lon = fibo_points[,1]))
mesh1 <- inla.mesh.2d(loc = fibo_points_xyz, cutoff = 0.01, max.edge = 0.5)
summary(mesh1)

## 5 Find GIA prior mean at mesh vertices
meshLL <- Lxyz2ll(list(x=mesh1$loc[,1], y = mesh1$loc[,2], z = mesh1$loc[,3]))
meshLL$lon <- ifelse(meshLL$lon >= -0.5, meshLL$lon,meshLL$lon + 360)
mesh_sp <- SpatialPoints(data.frame(lon = meshLL$lon, lat = meshLL$lat), proj4string = CRS("+proj=longlat")) 
mesh_idx <- over(mesh_sp, Plist)
GIA_prior <- ice6g$trend[mesh_idx]

## 6 Detrend the GPS data
GPS_data <- GPSV4b
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_data$lat, lon = GPS_data$lon))
GPS_sp <- SpatialPoints(data.frame(lon = ifelse(GPS_data$lon>359.5, GPS_data$lon - 360, GPS_data$lon), 
                                   lat = GPS_data$lat), proj4string = CRS("+proj=longlat"))
GPS_idx <- over(GPS_sp, Plist)
GPS_mu <- ice6g$trend[GPS_idx]
GPS_data$trend0 <- GPS_data$trend - GPS_mu


## 7 For the SPDE model prior used in the default setting of INLA, we need the following transformation 
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0
hyper <- list(prec = list(fixed = TRUE, initial = 0))
formula = y ~ -1 +  f(GIA, model = GIA_spde)

## 8 Load the pseudo polygon
zeroPolygon <- readOGR(dsn = "/./projects/GlobalMass/WP1-BHM/Experiment1b/shapefiles", layer = "zero03")
## Remove polygons that are too small
zeroPolys <- zeroPolygon@polygons[[1]]@Polygons
polyareas <- sapply(zeroPolys, function(x) x@area)
polyholes <- sapply(zeroPolys, function(x) x@hole)
zeropolys2 <- zeroPolys[polyareas > 200 ] 
zeroPoly <- zeroPolygon
zeroPoly@polygons[[1]]@Polygons <- zeropolys2

## 9 Generate mesh for parttion model
#### Dense points outside the polygons
pinPoly <- unlist(over(zeroPoly, SpatialPoints(coords = fibo_points), returnList=T))
fibo_inSub<- fibo_points[-pinPoly,]
#### Sparse points in the polygons
fibo_points2 <- fiboSphere(N = 500, L0=TRUE)
pinPoly <- unlist(over(zeroPoly, SpatialPoints(coords = fibo_points2), returnList=T))
fibo_inPoly<- fibo_points2[pinPoly,]
fibo_points_all <- rbind(fibo_inPoly, fibo_inSub)
mesh_points_xyz <- do.call(cbind, Lll2xyz(lat = fibo_points_all[,2], lon = fibo_points_all[,1]))
mesh2 <- inla.mesh.2d(loc = mesh_points_xyz, cutoff = 0.01, max.edge = 0.5)
mesh2 <- inla.mesh.2d(loc = mesh2$loc, cutoff = 0.01, max.edge = 0.5)

## 10 For this new mesh the GIA prior mean for the mesh vertices need to re-calculate.
meshLL <- Lxyz2ll(list(x=mesh2$loc[,1], y = mesh2$loc[,2], z = mesh2$loc[,3]))
meshLL$lon <- ifelse(meshLL$lon >= -0.5, meshLL$lon,meshLL$lon + 360)
mesh_sp <- SpatialPoints(data.frame(lon = meshLL$lon, lat = meshLL$lat), proj4string = CRS("+proj=longlat")) 
mesh_idx <- over(mesh_sp, Plist)
GIA_prior2 <- ice6g$trend[mesh_idx]

## 11 parttion the mesh2
mesh2 <- dt.mesh.addon.posTri(mesh = mesh2, globe = TRUE)
Tlonlat <- Lxyz2ll(list(x = mesh2$posTri[,1], y = mesh2$posTri[,2], z = mesh2$posTri[,3]))
Tlonlat$lon <- ifelse(Tlonlat$lon >=0, Tlonlat$lon, Tlonlat$lon + 359)
mesh2$Trill <- cbind(lon = Tlonlat$lon, lat =Tlonlat$lat)
TinPoly <- unlist(over(zeroPoly, SpatialPoints(coords=mesh2$Trill), returnList=T))
TAll <- 1:mesh2$t
ToutPoly <- TAll[-TinPoly]
Omega = dt.Omega(list(TinPoly, 1:mesh2$t), mesh2)


## 12 Prepare GPS data for non-stationary models 
## 12.1 Remove GPS data inside the polygon -- GPS data 2
GPS_inPoly <- unlist(over(zeroPoly, SpatialPoints(coords = cbind(GPSV4b$lon, GPSV4b$lat), 
                                                  proj4string = CRS(proj4string(zeroPoly))), returnList=T))
GPS_All <- 1:nrow(GPS_data)
GPS_outPoly <- GPS_All[-GPS_inPoly]
GPS_data2 <- GPS_data[GPS_outPoly,]
GPS_loc2 <- GPS_loc[GPS_outPoly,]

## 12.2 Pseudo GPS data along the boundary -- GPS data 3
boundlines <- as(zeroPoly, 'SpatialLines') 
obs_bounds <- spsample(boundlines, n = 30, type = "regular") # note points more than specified
proj4string(obs_bounds) <- proj4string(Plist)
pobs_idx <- over(obs_bounds, Plist)
GIA_pobs <- ice6g$trend[pobs_idx]
nobsb <-nrow(obs_bounds@coords) 
GPS_data3 <- data.frame(ID = rep("pseudo", nobsb), lon = obs_bounds@coords[,1], lat = obs_bounds@coords[,2],
                        trend = rep(0, nobsb), std = rep(0.05, nobsb), trend0 = -GIA_pobs)
GPS_loc3 <- do.call(cbind, Lll2xyz(lat = obs_bounds@coords[,2], lon = obs_bounds@coords[,1]))

## 12.3 Pseudo GPS data within the zero region -- sparse samples -- GPS data 4
## Create a sparse and regular GPS observation in side the zero polygon
proj4string(zeroPoly) <- "+proj=longlat"
GPS_ll <- spsample(zeroPoly, n = 30, type = "Fibonacci")
GIA_pobs <-  ice6g$trend[over(GPS_ll, Plist)]
nobsb <-nrow(GPS_ll@coords)
GPS_data4 <- data.frame(ID = rep("pseudo", nobsb), lon = GPS_ll@coords[,1], lat = GPS_ll@coords[,2],
                        trend = rep(0, nobsb), std = rep(0.05, nobsb), trend0 = -GIA_pobs)
GPS_loc4 <- do.call(cbind, Lll2xyz(lat = GPS_data4$lat, lon = GPS_data4$lon))


## 12.4 Pseudo GPS data within the zero region -- mesh vertices -- GPS data 5
## Find the mesh nodes incide the polygons
Vll <- Lxyz2ll(list(x = mesh2$loc[,1], y = mesh2$loc[,2], z = mesh2$loc[,3]))
Vll$lon <- ifelse(Vll$lon < 0, Vll$lon + 360, Vll$lon)
Vll <- cbind(Vll$lon, Vll$lat)
VinPoly <- unlist(over(zeroPoly, SpatialPoints(coords=Vll,proj4string = CRS(proj4string(zeroPoly))), returnList=T))
obs_inpoly <- Vll[VinPoly,]
obs_pseudo <-  SpatialPoints(coords = obs_inpoly)
proj4string(obs_pseudo) <- proj4string(Plist)
## Now assemble the new GPS data 
pobs_idx <- over(obs_pseudo, Plist)
GIA_pobs <- ice6g$trend[pobs_idx]
nobsb <-nrow(obs_pseudo@coords)
GPS_data5 <- data.frame(ID = rep("pseudo", nobsb), lon = obs_pseudo@coords[,1], lat = obs_pseudo@coords[,2],
                        trend = rep(0, nobsb), std = rep(0.05, nobsb), trend0 = -GIA_pobs)
GPS_loc5 <- mesh2$loc[VinPoly,]


##########################################################################
## 1 Stationary model
##########################################################################
## Build the SPDE model
mesh <- mesh1
gps_data <- GPS_data2
gps_loc <- GPS_loc2
GIA_spde <- inla.spde2.matern(mesh, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

## Link the process to observations and predictions
A_data <- inla.spde.make.A(mesh = mesh, loc = gps_loc)
A_pred <- inla.spde.make.A(mesh = mesh, loc = rbind(gps_loc, mesh$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=gps_data$trend0), A = list(A_data),
                     effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

## ----inla_run, message = FALSE, warning = FALSE--------------------------
## Fix the GPS errors
prec_scale <- c(1/gps_data$std^2, rep(1, nrow(A_pred)))

res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))

## ----inla_res, include = TRUE, cache=TRUE--------------------------------
INLA_pred <- res_inla$summary.linear.predictor
## Extract and project predictions
pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
GPS_idx <- pred_idx[1:nrow(gps_data)]
GIA_idx <- pred_idx[-(1:nrow(gps_data))]

## GPS 
GPS_u <- INLA_pred$sd[GPS_idx]
GPS_pred <- data.frame(lon = gps_data$lon, lat = gps_data$lat, u = GPS_u)

## GIA
GIA_diff <- INLA_pred$mean[GIA_idx] 
GIA_m <- GIA_diff + GIA_prior
GIA_u <- INLA_pred$sd[GIA_idx]
proj <- inla.mesh.projector(mesh, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-90, 90))
GIA_grid <- expand.grid(proj$x, proj$y)
GIA_predn <- data.frame(diff = GIA_diff, mean = GIA_m, u=GIA_u)
GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                       diff = as.vector(inla.mesh.project(proj, as.vector(GIA_diff))),
                       mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                       u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))),
                       model = rep("stationary", nrow(GIA_grid)))

ress1 <- list(res_inla = res_inla, spde = GIA_spde, st = stGIA, 
              mesh = mesh1, GPS_pred = GPS_pred, GIA_pred = GIA_pred, GIA_predn = GIA_predn)



##########################################################################
## 2 Subset model
##########################################################################
mesh_outPoly <- mesh.sub(mesh2, Omega, 2)
mesh <- mesh_outPoly
gps_data <- rbind(GPS_data2, GPS_data3)
gps_loc <- rbind(GPS_loc2, GPS_loc3)

GIA_spde <- inla.spde2.matern(mesh, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

## Link the process to observations and predictions
A_data <- inla.spde.make.A(mesh = mesh, loc =  rbind(gps_loc))
A_pred <- inla.spde.make.A(mesh = mesh, loc = rbind(gps_loc, mesh$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=c(gps_data$trend0)), A = list(A_data),
                     effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

## Fix the GPS errors
prec_scale <- c(1/gps_data$std^2, rep(1, nrow(A_pred)))

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))

## Extract and project predictions
INLA_pred <- res_inla$summary.linear.predictor
pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
GPS_idx <- pred_idx[1:nrow(gps_data)]
GIA_idx <- pred_idx[-c(1:nrow(gps_data))]

## GPS
GPS_u <- INLA_pred$sd[GPS_idx]
GPS_pred <- data.frame(lon = gps_data$lon, lat = gps_data$lat, u = GPS_u)

## GIA
GIA_diff <- INLA_pred$mean[GIA_idx]

## For this new mesh the GIA prior mean for the mesh vertices need to re-calculate.
meshLL <- Lxyz2ll(list(x=mesh_outPoly$loc[,1], y = mesh_outPoly$loc[,2], z = mesh_outPoly$loc[,3]))
meshLL$lon <- ifelse(meshLL$lon >= -0.5, meshLL$lon,meshLL$lon + 360)
mesh_sp <- SpatialPoints(data.frame(lon = meshLL$lon, lat = meshLL$lat), proj4string = CRS("+proj=longlat")) 
mesh_idx <- over(mesh_sp, Plist)
GIA_prior3 <- ice6g$trend[mesh_idx]

GIA_m <- GIA_diff + GIA_prior3
GIA_u <- INLA_pred$sd[GIA_idx]
proj <- inla.mesh.projector(mesh, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-90, 90))
GIA_grid <- expand.grid(proj$x, proj$y)
GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                       diff = as.vector(inla.mesh.project(proj, as.vector(GIA_diff))),
                       mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                       u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))),
                       model = rep("subset", nrow(GIA_grid)))
GIA_predn <- data.frame(diff = GIA_diff, mean = GIA_m, u=GIA_u)
ress2 <- list(res_inla = res_inla, spde = GIA_spde, st = stGIA,
              mesh = mesh, GPS_pred = GPS_pred, GIA_pred = GIA_pred, GIA_predn = GIA_predn)


##########################################################################
## 3 Constrained partition model
##########################################################################
mesh <- mesh2
gps_data <- rbind(GPS_data2, GPS_data4)
gps_loc <- rbind(GPS_loc2, GPS_loc4)

Q.mixture = dt.create.Q(mesh, Omega, fixed.ranges = c(4, NA))
prior0 <- list(sigma = tsigma, range = matrix(trho, ncol = 2))
log.prior <- dt.create.prior.log.norm(prior.param = prior0) 
GIA_spde = dt.inla.model(Q = Q.mixture, log.prior=log.prior)

## Link the process to observations and predictions
A_data <- inla.spde.make.A(mesh = mesh, loc =  gps_loc)
A_pred <- inla.spde.make.A(mesh = mesh, loc = rbind(gps_loc, mesh$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=gps_data$trend0), A = list(A_data),
                     effects = list(GIA = 1:mesh$n), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA = 1:mesh$n), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)
## Fix the GPS errors
prec_scale <- c(1/gps_data$std^2,  rep(1, nrow(A_pred)))

## Run INLA
res_inla <- inla(formula, data=inla.stack.data(stGIA), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A = inla.stack.A(stGIA), compute = TRUE))

## Extract and project predictions
INLA_pred <- res_inla$summary.linear.predictor
pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
GPS_idx <- pred_idx[1:nrow(gps_data)]
GIA_idx <- pred_idx[-c(1:nrow(gps_data))]

## GPS
GPS_u <- INLA_pred$sd[GPS_idx]
GPS_pred <- data.frame(lon = gps_data$lon, lat = gps_data$lat, u = GPS_u)

## GIA
GIA_diff <- INLA_pred$mean[GIA_idx]
GIA_m <- GIA_diff + GIA_prior2
GIA_u <- INLA_pred$sd[GIA_idx]
proj <- inla.mesh.projector(mesh, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-90, 90))
GIA_grid <- expand.grid(proj$x, proj$y)
GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                       diff = as.vector(inla.mesh.project(proj, as.vector(GIA_diff))),
                       mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                       u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))),
                       model = rep("mixture", nrow(GIA_grid)))
GIA_predn <- data.frame(diff = GIA_diff, mean = GIA_m, u=GIA_u)

ress3 <- list(res_inla = res_inla, spde = GIA_spde, st = stGIA,
              mesh = mesh, GPS_pred = GPS_pred, GIA_pred = GIA_pred, GIA_predn = GIA_predn)


save(ress1, ress2, ress3,  file ="/./projects/GlobalMass/WP1-BHM/Experiment1b/GIA_RGL/GIA_compare2.RData")
