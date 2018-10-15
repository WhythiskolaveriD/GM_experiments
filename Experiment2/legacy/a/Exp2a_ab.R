## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

wd <-"~/"
dd <-"/./projects/GlobalMass/"

library(rgdal); library(sp);library(GEOmap)
library(INLA)
library(ncdf4)
library(gstat)

source(paste0(wd,"gmrcode/BHM_sphere/functions.R"))
source(paste0(wd, "gmrcode/BHM_sphere/partition_fun.R"))

## Genereate Fibonacci points on the sphere
fibo_points <- fiboSphere(N = 12960, L0 = TRUE)
fibo_points_xyz <- do.call(cbind, Lll2xyz(lat = fibo_points[,2], lon = fibo_points[,1]))
mesh0 <- inla.mesh.2d(loc = fibo_points_xyz, cutoff = 0.01, max.edge = 1)
## Make this "smoother"
mesh0 <- inla.mesh.2d(loc = mesh0$loc, cutoff = 0.01, max.edge = 1)

## ----loaddata, include=FALSE, message=FALSE, warning=FALSE---------------
## Load Maike's GRACE data
gdata <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_trends_v03.txt"), header = T)
gm_loc <-  read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_loc_v03.txt"), skip = 1)
n_grace <- ncol(gm_loc)
n_ll <- (n_grace-4)/2
names(gm_loc) <- c("id", "area", "lon_c", "lat_c", paste0(c("lon", "lat"), rep(1:n_ll, each = 2)))

## Create spatial polygons data frame
Polygon_list <- list()
for(i in 1:nrow(gm_loc)){
  lons <- na.omit(as.numeric(c(gm_loc[i, seq(5,ncol(gm_loc), 2)], gm_loc[i, 5])))
  lats <- na.omit(as.numeric(c(gm_loc[i, seq(6,ncol(gm_loc), 2)], gm_loc[i, 6])))
  Polygon_list[[i]] <- Polygon(cbind(lons, lats))
}

Polygons_list <- lapply(1:length(Polygon_list), function(x) Polygons(list(Polygon_list[[x]]), x))
SpPolygon <- SpatialPolygons(Polygons_list, proj4string = CRS("+proj=longlat"))

graceM <- SpatialPolygonsDataFrame(SpPolygon,gdata)
graceM$area <- gm_loc$area
gm_area <- geosphere::areaPolygon(graceM)/(1000^2) # km^2
graceM$area2 <- gm_area
graceM$diff <- graceM$area/graceM$area2
graceM$lon_c <- gm_loc$lon_c
graceM$lat_c <- gm_loc$lat_c

## Load the ice6g-ewh
ice6g <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5_ewh.txt"), header = T)
names(ice6g)[2:3] <- c("lon", "lat")
ice6g$lon <- ifelse(ice6g$lon <0, ice6g$lon + 360, ice6g$lon)
ice6g<- ice6g[order(ice6g$lat,ice6g$lon ),] # re-order the data according to the coordinates
coordinates(ice6g) <-c("lon", "lat")
gridded(ice6g) <- TRUE
proj4string(ice6g) <- CRS("+proj=longlat")
ice6g$areas <- geosphere::areaPolygon(as(ice6g, "SpatialPolygons"))/(1000^2)



## Load Rory's GRACE data
gr1 <- read.table(paste0(dd, "WP1-BHM/Experiment2a/grace_grid.txt"))
gr2 <-read.table(paste0(dd, "WP1-BHM/Experiment2a/grace-gia_grid.txt"))
lon <- 0.5:359.5
lat <- -89.5:89.5
locs <- expand.grid(lon, lat)
graceR <- data.frame(lon = locs[,1], lat = locs[,2], trend = gr1$V1, trendMgia = gr2$V1)
coordinates(graceR) <-c("lon", "lat")
gridded(graceR) <- TRUE
proj4string(graceR) <- CRS("+proj=longlat")
graceR$areas <- geosphere::areaPolygon(as(graceR, "SpatialPolygons"))/(1000^2)

## Add GRACE - i6g trend by hand
graceR$trendMgia2 <- graceR$trend - ice6g$trend
graceR$giadiff <- graceR$trendMgia - graceR$trendMgia2


## Load Ocean boundary
Ocean <- readOGR(dsn = paste0(dd, "WP1-BHM/maps/ne_110m_ocean"), layer = "ne_110m_ocean")

## ----gerror, echo=FALSE--------------------------------------------------
graceR$std <- over(graceR, graceM)$std

## ----Grace_prior, echo = FALSE-------------------------------------------
#v2 <- variogram(trendMgia~1, graceR[sample(1:nrow(graceR), 3000),]) 
#plot(v2)

## Priors mean and variance for the parameters: rho and sigma
mu_r <- 2500/6371
v_r <- 1
mu_s <- 25 # area scale to about 100km^2 ~ 1 degree resolution
v_s <- 50^2

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


## ----grace_area, cache = TRUE--------------------------------------------
grace_loc <- do.call(cbind, Lll2xyz(lon = graceR@coords[,1], lat = graceR@coords[,2]))
A_GRACE_data <- inla.spde.make.A(mesh = mesh0, loc = grace_loc)
A_M_pred <- inla.spde.make.A(mesh = mesh0, loc = rbind(mesh0$loc))

## ----stacks, eval = FALSE------------------------------------------------
## ## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=graceR$trendMgia), A = list(A_GRACE_data),
                     effects = list(M= 1:M_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(rbind(A_GRACE_data, A_M_pred)),
                      effects = list(M=1:M_spde$n.spde), tag = "pred")
stM <- inla.stack(st.est, st.pred)

## ----inla_run_grace, include = TRUE, eval = FALSE------------------------
## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))
prec_scale <- c(1/graceR$std^2, rep(1, nrow(A_GRACE_data) + nrow(A_M_pred)))


## Add constraints sum to zero
## The formular for modelling the SSH mean
formula = y ~ -1 +  f(M, model = M_spde, extraconstr = list(A = matrix(1, nrow = 1, ncol = mesh0$n), e = 0))

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stM, spde = M_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stM), compute =TRUE))

## ----inla_res_grace, include = TRUE, eval = FALSE, echo = FALSE----------
INLA_pred <- res_inla$summary.linear.predictor
## Extract and project predictions
pred_idx <- inla.stack.index(stM, tag = "pred")$data[1:64800]


## mass
M_m <- INLA_pred$mean[pred_idx]
M_u <- INLA_pred$sd[pred_idx]
proj <- inla.mesh.projector(mesh0, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-89.5, 89.5))
M_grid <- expand.grid(proj$x, proj$y)
M_pred <- data.frame(lon = M_grid[,1], lat = M_grid[,2],
                     mean = M_m, u = M_u, diff = M_m - graceR@data$trendMgia)

res_M2 <- list(res_inla = res_inla, spde = M_spde, st = stM,
            mesh = mesh0,  M_pred = M_pred, Adata = A_GRACE_data, Apred = A_M_pred)

graceR@data$pred_mean <- M_m
graceR@data$pred_u <- M_u
grace@data$pred_diff <- M_m - graceR@data$trendMgia

save(res_M2, graceR, file =paste0(dd, "WP1-BHM/Experiment2a/exp2a_M2zero.RData"))
