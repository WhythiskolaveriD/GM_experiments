## ----ini_poly, message = FALSE, warning=FALSE----------------------------
## load library and functions
library(sp); library(INLA); library(GEOmap); library(rgdal)
library(ggplot2); library(grid); library(gridExtra)
if(Sys.info()['sysname'] == "Linux"){INLA:::inla.dynload.workaround()}
source("functions.R")
source("functions-barriers-dt-models-march2017.R")

## Load the pseudo polygon
#### 1 Load GIA prior
if(Sys.info()["sysname"] == "Windows"){
  zeroPolygon <- readOGR(dsn = "Z:/WP1-BHM/Experiment1b/shapefiles", layer = "zero03")
}else if(Sys.info()["sysname"] == "Ubuntu"){
  zeroPolygon <- readOGR(dsn = "Z:/WP1-BHM/Experiment1b/shapefiles", layer = "zero03")
}else{
  zeroPolygon <- readOGR(dsn = "/./projects/GlobalMass/WP1-BHM/Experiment1b/shapefiles", layer = "zero03")
}
## Remove polygons that are too small
zeroPolys <- zeroPolygon@polygons[[1]]@Polygons
polyareas <- sapply(zeroPolys, function(x) x@area)
polyholes <- sapply(zeroPolys, function(x) x@hole)

zeropolys2 <- zeroPolys[polyareas > 200 ] 
zeroPoly <- zeroPolygon
zeroPoly@polygons[[1]]@Polygons <- zeropolys2
plot(zeroPoly, col = "blue", main = "zero regions -- threshold = 0.3 (refined)")

## ----mesh, include=TRUE, cache=TRUE--------------------------------------
#### Dense points over in the subset
fibo_points <- fiboSphere(N = 12960, L0 = TRUE)
pinPoly <- unlist(over(zeroPoly, SpatialPoints(coords = fibo_points), returnList=T))
fibo_inSub<- fibo_points[-pinPoly,]
plot(zeroPoly)
points(fibo_inSub, pch = ".")

#### Sparse points in the polygons
fibo_points <- fiboSphere(N = 500, L0=TRUE)
pinPoly <- unlist(over(zeroPoly, SpatialPoints(coords = fibo_points), returnList=T))
fibo_inPoly<- fibo_points[pinPoly,]
plot(zeroPoly)
points(fibo_inPoly, pch = ".")
points(fibo_inSub, pch = ".")

fibo_points_all <- rbind(fibo_inPoly, fibo_inSub)
library(GEOmap);library(INLA)
if(Sys.info()['sysname'] == "Linux"){INLA:::inla.dynload.workaround()}
mesh_points_xyz <- do.call(cbind, Lll2xyz(lat = fibo_points_all[,2], lon = fibo_points_all[,1]))
mesh <- inla.mesh.2d(loc = mesh_points_xyz, cutoff = 0.01, max.edge = 0.5)
summary(mesh) # give the desired number of vertices and triangles.

## ----mesh2, include = TRUE, cache=TRUE-----------------------------------
mesh <- dt.mesh.addon.posTri(mesh = mesh, globe = TRUE)
Tlonlat <- Lxyz2ll(list(x = mesh$posTri[,1], y = mesh$posTri[,2], z = mesh$posTri[,3]))
Tlonlat$lon <- ifelse(Tlonlat$lon >=0, Tlonlat$lon, Tlonlat$lon + 359)
mesh$Trill <- cbind(lon = Tlonlat$lon, lat =Tlonlat$lat)
TinPoly <- unlist(over(zeroPoly, SpatialPoints(coords=mesh$Trill), returnList=T))
TAll <- 1:mesh$t
ToutPoly <- TAll[-TinPoly]
Omega = dt.Omega(list(TinPoly, 1:mesh$t), mesh)

## Plot the result in 3d
plot(mesh, t.sub = Omega[[2]], col = "lightblue", rgl = TRUE )
plot(mesh, t.sub = Omega[[1]], col = "yellow",  rgl = TRUE, add = TRUE)

## ----mesh4, include=TRUE, cache=TRUE-------------------------------------
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0


GIA_spde <- inla.spde2.matern(mesh, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

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

meshLL <- Lxyz2ll(list(x=mesh$loc[,1], y = mesh$loc[,2], z = mesh$loc[,3]))
meshLL$lon <- ifelse(meshLL$lon >= -0.5, meshLL$lon,meshLL$lon + 360)
mesh_sp <- SpatialPoints(data.frame(lon = meshLL$lon, lat = meshLL$lat), proj4string = CRS("+proj=longlat")) 
mesh_idx <- over(mesh_sp, Plist)
GIA_prior <- ice6g$trend[mesh_idx]

#### 2 Load GPS data
if(Sys.info()["sysname"] == "Windows"){
  GPSV4b <- read.table("Z:/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)
}else if(Sys.info()["sysname"] == "Ubuntu"){
  GPSV4b <- read.table("~/globalmass/WP2-SolidEarth/GPS/NGL/BHMinputFiles/GPS_v04b.txt", header = T)
}else{
  GPSV4b <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)
}

## ----data, include = TRUE, cache=TRUE------------------------------------
GPS_inPoly <- unlist(over(zeroPoly, SpatialPoints(coords = cbind(GPSV4b$lon, GPSV4b$lat)), returnList=T))
GPS_All <- 1:nrow(GPSV4b)
GPS_outPoly <- GPS_All[-GPS_inPoly]
plot(GPSV4b[GPS_outPoly,c("lon", "lat")], pch = "+")

GPS_data <- GPSV4b[GPS_outPoly,]
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_data$lat, lon = GPS_data$lon))
GPS_sp <- SpatialPoints(data.frame(lon = ifelse(GPS_data$lon>359.5, GPS_data$lon - 360, GPS_data$lon), 
                                   lat = GPS_data$lat), proj4string = CRS("+proj=longlat"))

GPS_idx <- over(GPS_sp, Plist)
GPS_mu <- ice6g$trend[GPS_idx]
GPS_data$trend0 <- GPS_data$trend - GPS_mu

## ----data2, include = TRUE, cache=TRUE-----------------------------------
## get the boundary of the polygons
boundlines <- as(zeroPoly, 'SpatialLines') 
obs_bounds <- spsample(boundlines, n = 50, type = "regular") # note points more than specified
## The mesh nodes incide the polygons
obs_inpoly <- SpatialPoints(coords = mesh$Trill[TinPoly,])
obs_pseudo <- rbind(obs_bounds, obs_inpoly)
proj4string(obs_pseudo) <- proj4string(Plist)

## Find the ice6g values
pobs_idx <- over(obs_pseudo, Plist)
GIA_pobs <- ice6g$trend[pobs_idx]
nb <- length(obs_bounds)
ni <- length(obs_inpoly)
## the average of these observations 
GIA_Poly_obs <- -mean(GIA_pobs[-(1:nb)])


nobsb <-nrow(obs_pseudo@coords)
obs_df <- data.frame(ID = rep("pseudo", nobsb), lon = obs_pseudo@coords[,1], lat = obs_pseudo@coords[,2],
                     trend = rep(0, nobsb), std = rep(0.05,nobsb), trend0 = -GIA_pobs)
obs_xyz <- do.call(cbind, Lll2xyz(lat = obs_pseudo@coords[,2], lon = obs_pseudo@coords[,1]))

GPS_all <- rbind(GPS_data, obs_df)
GPS_all_loc <- rbind(GPS_loc, obs_xyz)

## ----inla, include = TRUE, cache = TRUE----------------------------------
## Link the process to observations and predictions
Mesh_GIA <- mesh
npts <- nrow(GPS_data) + nb
block <- c(1:npts, rep(npts+1, ni))
A_data <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_all_loc, 
                           block = block, n.block = npts + 1, block.rescale = "count")
A_pred <- inla.spde.make.A(mesh = Mesh_GIA, loc = rbind(GPS_all_loc, Mesh_GIA$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=c(GPS_data$trend0, obs_df$trend0[1:nb], GIA_Poly_obs)), A = list(A_data),
                     effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

## Fix the GPS errors
hyper <- list(prec = list(fixed = TRUE, initial = 0))
formula = y ~ -1 +  f(GIA, model = GIA_spde)
prec_scale <- c(1/GPS_data$std^2, rep(1/0.05^2, nb), 1/0.01^2, rep(1, nrow(A_pred)))

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))

save(res_inla, file = "/./projects/GlobalMass/WP1-BHM/Experiment1b/GIA_RGL/res3.RData")


## ----inla_res, include = TRUE, cache=TRUE--------------------------------
INLA_pred <- res_inla$summary.linear.predictor

## Extract and project predictions
pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
GPS_idx <- pred_idx[1:(nrow(GPS_data)+1)]
GIA_idx <- pred_idx[-c(1:(nrow(GPS_all)+1))]

## GPS 
GPS_u <- INLA_pred$sd[GPS_idx] 
GPS_pred <- data.frame(lon = GPS_all$lon, lat = GPS_all$lat, u = GPS_u)

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
                    colpal = colorRamps::matlab.like(12),  title = "Predicted uncertainties")

## Display 
print(map_prior)
print(map_GIA)
print(map_diff)
print(map_sd)

