library(rgdal); library(sp);library(GEOmap)
library(INLA)
source("C:/ZSwork/glbm/BHM_sphere/functions.R")
source("C:/ZSwork/glbm/BHM_sphere/partition_fun.R")
## Genereate Fibonacci points on the sphere
fibo_points <- fiboSphere(N = 12960, L0 = TRUE)
fibo_points_xyz <- do.call(cbind, Lll2xyz(lat = fibo_points[,2], lon = fibo_points[,1]))
mesh0 <- inla.mesh.2d(loc = fibo_points_xyz, cutoff = 0.01, max.edge = 1)
## Make this "smoother"
mesh0 <- inla.mesh.2d(loc = mesh0$loc, cutoff = 0.01, max.edge = 1)

## ----grace_data----------------------------------------------------------
grace_data <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_trends_v03.txt", header = T)
grace_loc <-  read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_loc_v03.txt", skip = 1)
n_grace <- ncol(grace_loc)
n_ll <- (n_grace-4)/2
names(grace_loc) <- c("id", "area", "lon_c", "lat_c", paste0(c("lon", "lat"), rep(1:n_ll, each = 2)))

## Create spatial polygons data frame
Polygon_list <- list()
for(i in 1:nrow(grace_loc)){
  lons <- na.omit(as.numeric(c(grace_loc[i, seq(5,ncol(grace_loc), 2)], grace_loc[i, 5])))
  lats <- na.omit(as.numeric(c(grace_loc[i, seq(6,ncol(grace_loc), 2)], grace_loc[i, 6])))
  Polygon_list[[i]] <- Polygon(cbind(lons, lats))
}

Polygons_list <- lapply(1:length(Polygon_list), function(x) Polygons(list(Polygon_list[[x]]), x))
SpPolygon <- SpatialPolygons(Polygons_list, proj4string = CRS("+proj=longlat"))

grace_sp <- SpatialPolygonsDataFrame(SpPolygon,grace_data)

## Priors mean and variance for the parameters: rho and sigma
mu_r <- 2000/6371
v_r <- 1
mu_s <- 2
v_s <- 5^2

## Transform the parameters for the SPDE_GMRF approximation
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

## Build the SPDE model with the prior
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

M_spde <- inla.spde2.matern(mesh0, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                               theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))


## ----grace_link----------------------------------------------------------
## For each polygon observation we generate the regular spaced grid and the number of grid cell is proportional to the area of the polygon
grace_area <- geosphere::areaPolygon(grace_sp)/(1000^2)
grace_sp$area <-grace_area
area_mean <- mean(grace_area)

## Generate the integration grid for each polygons
poly_block <- function(i, dis = 10){
  sp_i <- SpatialPolygons(list(grace_sp@polygons[[i]]), proj4string=CRS("+proj=longlat"))
  area_i <- grace_sp$area[i]
  grid_i <- spsample(sp_i, n = round(area_i/dis^2), type = "regular", offset=c(0.5, 0.5))
  ngrid_i <- length(grid_i)
  grid_xyz <- do.call(cbind, Lll2xyz(lat = grid_i@coords[,2], lon = grid_i@coords[,1]))
  block_i <- rep(i, ngrid_i) 
  weights <- rep(area_i/ngrid_i, ngrid_i)
  return(list(grid_xyz = grid_xyz, block = block_i, weights = weights, ngrid = ngrid_i))
}

grace_block <- lapply(1:nrow(grace_sp), poly_block, dis = 10)

grid_xyz <- do.call(rbind, lapply(grace_block, "[[", "grid_xyz"))
grid_block <- do.call(c, lapply(grace_block, "[[", "block"))
weights <- do.call(c, lapply(grace_block, "[[", "weights"))

A_GRACE_data <- inla.spde.make.A(mesh = mesh0, loc = grid_xyz, block = grid_block,  weights = weights)

## Same for prediction on a 1 degree resolution grid, we need to know the area of the grid for integration
## generate the prediction grid
gx <- seq(0, 359, 1)
gy <- seq(-89.5, 89.5)
grid_ll <- expand.grid(gx, gy)
pred_data <- data.frame(lon = grid_ll[,1], lat = grid_ll[,2])
coordinates(pred_data) <-c("lon", "lat")
gridded(pred_data) <- TRUE
pred_data <- as(pred_data, "SpatialPolygons")
proj4string(pred_data) <- CRS("+proj=longlat")
areas <- geosphere::areaPolygon(pred_data)/(1000^2)
grid_pred <- do.call(cbind,Lll2xyz(lat = grid_ll[,2], lon = grid_ll[,1]))
A_M_pred <- inla.spde.make.A(mesh = mesh0, loc = grid_pred)

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=grace_sp$mmweq * grace_sp$area), A = list(A_GRACE_data),
                     effects = list(M = 1:M_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(rbind(A_GRACE_data, A_M_pred)),
                      effects = list(M=1:M_spde$n.spde), tag = "pred")
stM <- inla.stack(st.est, st.pred)

## ----inla_run_grace, include = TRUE, eval = FALSE------------------------
## ## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))
prec_scale <- c(1/grace_sp$std^2, rep(1, nrow(A_M_pred) + nrow(A_GRACE_data)))

## The formular for modelling the SSH mean
formula = y ~ -1 +  f(M, model = M_spde)

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stM, spde = M_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stM), compute =TRUE))

## ----inla_res_grace, include = TRUE, eval = FALSE------------------------
INLA_pred <- res_inla$summary.linear.predictor
## Extract and project predictions
pred_idx <- inla.stack.index(stM, tag = "pred")$data
idx_grace <- pred_idx[1:nrow(A_GRACE_data)]
idx_grid <- pred_idx[-(1:nrow(A_GRACE_data))]

## SSH
M_m <- INLA_pred$mean[idx_grid]
M_u <- INLA_pred$sd[idx_grid]
proj <- inla.mesh.projector(mesh0, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-89.5, 89.5))
M_grid <- expand.grid(proj$x, proj$y)
M_pred <- data.frame(lon = M_grid[,1], lat = M_grid[,2],
                        mean = M_m,
                        u = M_u)

res_M <- list(res_inla = res_inla, spde = M_spde, st = stM, 
                 mesh = mesh0,  M_pred = M_pred,  Adata = A_GRACE_data, Apred = A_M_pred)

grace_m <- INLA_pred$mean[idx_grace]/grace_sp$area
grace_u <- INLA_pred$sd[idx_grace]
grace_sp$pred_mean <- grace_m
grace_sp$pred_u <- grace_u


save(res_M, grace_sp , file ="/./projects/GlobalMass/WP1-BHM/Experiment2a/exp2a_M.RData")
