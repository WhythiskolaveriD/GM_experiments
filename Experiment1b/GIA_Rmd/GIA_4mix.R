## ----load, message = FALSE-----------------------------------------------
## load library and functions
library(INLA)
library(sp); library(GEOmap); library(rgdal)
library(ggplot2); library(grid); library(gridExtra)
source("functions.R")
source("functions-barriers-dt-models-march2017.R")

## ----polygons, include = FALSE, message = FALSE, cache = TRUE------------
## Load the pseudo polygon
#### 1 Load GIA prior
zeroPolygon <- readOGR(dsn = "/./projects/GlobalMass/WP1-BHM/Experiment1b/shapefiles", layer = "zero03")

## Remove polygons that are too small
zeroPolys <- zeroPolygon@polygons[[1]]@Polygons
polyareas <- sapply(zeroPolys, function(x) x@area)
polyholes <- sapply(zeroPolys, function(x) x@hole)

zeropolys2 <- zeroPolys[polyareas > 200 ] 
zeroPoly <- zeroPolygon
zeroPoly@polygons[[1]]@Polygons <- zeropolys2

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
plot(mesh, t.sub = Omega[[2]])
plot(mesh, t.sub = Omega[[1]])

## ----mesh4, include=TRUE, cache=TRUE-------------------------------------
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

## Transform the parameters for the SPDE_GMRF approximation
Tlognorm <- function(mu, v){
  logv <- log(1 + v/mu^2)
  logmu <- log(mu^2) - 0.5*log(mu^2 + v)
  return(c(logmu, logv))
}
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0


## Create the precsion matrix function
Q.function = dt.create.Q(mesh, Omega)
ranges = c(log(0.001), lrho0)
# - the first range is for the barrier area
# - - it is not sensitive to the exact value here,
# - the second range is for the normal area
Q = Q.function(theta = c(lsigma0, ranges))
# - the precision matrix for fixed ranges

## Simulate the field
u = inla.qsample(n=1, Q=Q, seed = 2017)
u = u[ ,1]

## Plot the simulated field
local.plot.field(u, mesh = mesh, main="The true (simulated) spatial field", asp = 1)


## ----load_data0, include=FALSE, eval = TRUE, cache = TRUE----------------
#### 1 Load GIA prior
ice6g <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)


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
GPSV4b <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)


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

nobsb <-nrow(obs_pseudo@coords)
obs_df <- data.frame(ID = rep("pseudo", nobsb), lon = obs_pseudo@coords[,1], lat = obs_pseudo@coords[,2],
                     trend = rep(0, nobsb), std = rep(0.05,nobsb), trend0 = -GIA_pobs)
obs_xyz <- do.call(cbind, Lll2xyz(lat = obs_pseudo@coords[,2], lon = obs_pseudo@coords[,1]))

GPS_all <- rbind(GPS_data, obs_df)
GPS_all_loc <- rbind(GPS_loc, obs_xyz)

## ----inla, include = TRUE, cache = TRUE----------------------------------
Mesh_GIA <- mesh
A_data <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_all_loc)
A_pred <- inla.spde.make.A(mesh = Mesh_GIA, loc = rbind(GPS_all_loc, Mesh_GIA$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=GPS_all$trend0), A = list(A_data),
                     effects = list(GIA = 1:Mesh_GIA$n), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:Mesh_GIA$n), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

## Fix the GPS errors
hyper <- list(prec = list(fixed = TRUE, initial = 0))
prec_scale <- c(1/GPS_all$std^2, rep(1, nrow(A_pred)))


Q.mixture = dt.create.Q(mesh, Omega, fixed.ranges = c(0.001, NA))
# - We fix the barrier range to a different value than we 
#   used for simulations
# - - Why? It does not matter, as long as it is 'small' 
#     the models are very
#     similar
# - - This shows that you do not need to know the 
#     true 'barrier range'!

log.prior = dt.create.prior.log.exp(prior.param = c(1,1))
# - The prior parameters are the lambdas in the exponential 
#   priors for standard deviation and inverse-range

GIA_spde = dt.inla.model(Q = Q.mixture, log.prior=log.prior)

formula <- y ~ -1 + f(GIA, model=GIA_spde)
# - The spatial model component is different from before
# - The rest of the model setup is the same! 
#   (as in the stationary case)
# - - e.g. the inla(...) call below is the same, 
#     only this formula is different

## ----run_inla, include = TRUE, eval = FALSE------------------------------
res_inla <- inla(formula, data=inla.stack.data(stGIA), family = "gaussian",
                   scale =prec_scale, control.family = list(hyper = hyper),
                   control.predictor=list(A = inla.stack.A(stGIA), compute = TRUE))
save(res_inla, file = "/./projects/GlobalMass/WP1-BHM/Experiment1b/GIA_RGL/res4.RData")