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

grace_bdv <- data.frame(trend = as.numeric(trend_grace), std = as.numeric(err_grace), lon= lon, lat = lat)
coordinates(grace_bdv) <- c("lon", "lat")

## Load BDV's ICE6G-VM5 in ewh mm/yr
gia_ewh <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA_BDV/GIA_ewh_ICE6G"), sep = ",", head = TRUE)
gia_ewh$lat <- -gia_ewh$lat
names(gia_ewh) <- c("lat", "lon", "trend")
gia_ewh$trend <- gia_ewh$trend*1000


## Load the Ocean polygon (-180, 180)
Ocean <- readOGR(dsn = paste0(dd, "WP1-BHM/maps/ne_110m_ocean"), layer = "ne_110m_ocean")

## The ICE6G-VM5 in vlm mm/yr
ice6g_vlm <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt"), header = T)
names(ice6g_vlm)[2:3] <- c("lon", "lat")
ice6g_vlm$lon <- ifelse(ice6g_vlm$lon >180, ice6g_vlm$lon - 360, ice6g_vlm$lon)
coordinates(ice6g_vlm) <-c("lon", "lat")
gridded(ice6g_vlm) <- TRUE
proj4string(ice6g_vlm) <- CRS(proj4string(Ocean))
ice6g_vlm$areas <- geosphere::areaPolygon(as(ice6g_vlm, "SpatialPolygons"))/(1000^2) # calculated the grid areas
Oid <- unlist(over(Ocean, as(ice6g_vlm, "SpatialPoints"), returnList = T)) # Find the grid in Ocean
vlm <- sum(ice6g_vlm$trend[Oid]*ice6g_vlm$areas[Oid]/sum(ice6g_vlm$areas[Oid])) # the average vlm in the ocean
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
grace_bdv$trendgia <- grace_bdv$trend - gia_sp$trend[giaG]
mean(grace_bdv$trendgia)

## ----variogram, message = FALSE, cache = TRUE, echo = FALSE--------------
## learn mass parameters from GRACE data
#v_grace <- gstat::variogram(trendgia ~ 1, data = grace_bdv)

#plot(v_grace)
# -- the range is likely to be about 1500km 
# -- the variance is likely to be around 15mm/yr

# ## learn steric parameters from altimetry - grace - vlm
# ## sample grace value for ssh
# steric_sp <- SpatialPoints(coords = cbind(alt_data$lon, alt_data$lat), proj4string = CRS(proj4string(G_sp)))
# G_trend <- over(steric_sp, G_sp)$trendgia
# steric_learn <- na.omit(data.frame(lon = alt_data$lon, lat = alt_data$lat, trend = alt_data$trend_ssh - G_trend - vlm))
# 
# coordinates(steric_learn) <- c("lon", "lat")
# proj4string(steric_learn) <- CRS("+proj=longlat")
# v_steric <- lapply(rep(5e3, 10), function(x) {gstat::variogram( trend ~ 1, data = steric_learn[sample(1:nrow(steric_learn), x), ])})
# plot(v_steric[[1]]$gamma ~ v_steric[[1]]$dist, type = "l", xlab = "distance", ylab = "semivariance", 
#      ylim = c(0, 350), xlim = c(0, 6000))
# for(i in 2:10){
#   lines(v_steric[[i]]$gamma ~ v_steric[[i]]$dist , col = i)
# }
# -- the range is likely to be 2000km
# -- the variance is likely to be 15mm/yr

## ----prior---------------------------------------------------------------
## Priors for mass
mu_r1 <- 1500/6371 # mean of rho
v_r1 <- 1 # vague variace for rho
mu_s1 <- 15 # mean of sigma
v_s1 <- 30^2 # vaue variance for sigma

## Priors for steric
mu_r2 <- 2000/6371 # mean of rho
v_r2 <- 1 # vague variace for rho
mu_s2 <- 15 # mean of sigma
v_s2 <- 30^2 # vaue variance for sigma

## Transform the parameters for a lognormal distribution
## -- mass
trho1 <- Tlognorm(mu_r1, v_r1)
tsigma1 <- Tlognorm(mu_s1, v_s1)

lsigma1 <- tsigma1[1]
theta1_ss <- tsigma1[2]
lrho1 <- trho1[1]
theta1_rs <- trho1[2]
lkappa1 <- log(8)/2 - lrho1
ltau1 <- 0.5*log(1/(4*pi)) - lsigma1 - lkappa1

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
TinOcean <- unlist(over(Ocean, SpatialPoints(coords=mesh0$Trill, proj4string = CRS(proj4string(Ocean))), returnList=T))
TAll <- 1:mesh0$t
ToutOcean <- TAll[-TinOcean]
Omega = dt.Omega(list(TinOcean, 1:mesh0$t), mesh0)

mesh_ocean <- mesh.sub(mesh0, Omega, 1)
summary(mesh_ocean)

## ----spde----------------------------------------------------------------
mass_spde <- inla.spde2.matern(mesh0, B.tau = matrix(c(ltau1, -1, 1),1,3), B.kappa = matrix(c(lkappa1, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_ss), sqrt(1/theta1_rs)))

steric_spde <- inla.spde2.matern(mesh_ocean, B.tau = matrix(c(ltau2, -1, 1),1,3), B.kappa = matrix(c(lkappa2, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta2_ss), sqrt(1/theta2_rs)))


## ----link_data-----------------------------------------------------------
## The data
ygrace <- grace_bdv$trendgia
yaltv<- alt_data2$trend_vlm

## The A matrice for linking observation and processes

## Link GRACE
grace_loc <- do.call(cbind, Lll2xyz(lat = grace_bdv@coords[,2], lon = grace_bdv@coords[,1]))
A_grace <- inla.spde.make.A(mesh = mesh0, loc = grace_loc)

## Link Altimetry
A_altv1 <- inla.spde.make.A(mesh = mesh0, loc = alt_loc)
A_altv2 <- inla.spde.make.A(mesh = mesh_ocean, loc = alt_loc)

## The A matrices for predict ssh
A_grace2ocean <- inla.spde.make.A(mesh = mesh0, loc = mesh_ocean$loc)

## The errors
prec_scale <- c(1/grace_bdv$std^2, 1/alt_data2$err_ssh^2, rep(1, nrow(A_grace2ocean)))


## Create the mass and steric stacks
stkmass <- inla.stack(data = list(y=ygrace), A = list(A_grace),
                     effects = list(mass = 1:mass_spde$n.spde), tag = "mass")

stksteric <- inla.stack(data = list(y=yaltv), A = list(A_altv1, A_altv2),
                     effects = list(list(mass = 1:mass_spde$n.spde),
                                    list(steric = 1:steric_spde$n.spde)), tag = "steric")

## Create the stack for predicting SSH
stkssh <- inla.stack(data = list(y=NA), A = list(A_grace2ocean, 1),
                     effects = list(list(mass = 1:mass_spde$n.spde),
                                    list(steric = 1:steric_spde$n.spde)), tag = "ssh")
## join the two stacks
stkall <- inla.stack(stkmass, stksteric, stkssh)

## ----inla_run, include = TRUE, eval = FALSE------------------------------
## ## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))

## The formular -- we add the constraint that mass change sum to zero
formula = y ~ -1 + f(mass, model = mass_spde, extraconstr = list(A = matrix(1, nrow = 1, ncol = mesh0$n), e = 0)) +
  f(steric, model = steric_spde)

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stkall), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stkall), compute =TRUE))

## ----inla_res, include = TRUE, eval = FALSE, echo = FALSE----------------
## idices of the predicted processes
mass_idx <- inla.stack.index(stkall, tag = "mass")$effects
mass_idxg <- inla.stack.index(stkall, tag = "mass")$data
steric_idx <- inla.stack.index(stkall, tag = "steric")$effects[-c(1:mass_spde$n.spde)]
ssh_idx <- inla.stack.index(stkall, tag = "ssh")$data


## Extract the predictions
INLA_mass <- res_inla$summary.random$mass
INLA_steric <- res_inla$summary.random$steric
INLA_ssh <- res_inla$summary.linear.predictor[ssh_idx,]
INLA_grace <- res_inla$summary.linear.predictor[mass_idxg,]

## mass
mass_mean <- INLA_mass$mean
mass_u <- INLA_mass$sd
proj <- inla.mesh.projector(mesh0, projection = "longlat", dims = c(360,180), xlim = c(0.5,359.5), ylim = c(-89.5, 89.5))
mass_grid <- expand.grid(proj$x, proj$y)
mass_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                       mean = as.vector(inla.mesh.project(proj, as.vector(mass_mean))),
                       u = as.vector(inla.mesh.project(proj, as.vector(mass_u))))

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

ress_2b_checkb <- list(res_inla = res_inla, st = stkall,
                spde = list(mass = mass_spde, steric = steric_spde),
                mesh = list(mass = mesh0, steric = mesh_ocean),
                pred = list(mass = mass_pred, steric = steric_pred, ssh = ssh_pred, grace_pred = grace_pred))

save(ress_2b_checkb, file = "/./projects/GlobalMass/WP1-BHM/Experiment2b/Exp2b_checkb.RData")
