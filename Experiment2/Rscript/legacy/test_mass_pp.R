### Test script for partition mass


## Setup and load data
############################################################################################################
wd <- "~/"
dd <- "/./projects/GlobalMass/"

library(rgdal); library(sp);library(GEOmap)
library(INLA)
library(ncdf4)
source(paste0(wd,"gmrcode/BHM_sphere/functions.R"))
source(paste0(wd, "gmrcode/BHM_sphere/partition_fun.R"))
source(paste0(wd, "gmrcode/BHM_sphere/pp_functions.R"))


## Load BDV's GRACE data -- note the data are not evenly spaced!
grace_nc <- nc_open(paste0(dd, "WP2-SolidEarth/Bramha/BHMinput/GRACE_trends_Oceanadded_Bramha.nc"))
print(grace_nc)
lon <- ncvar_get(grace_nc, "Long")
lat <- ncvar_get(grace_nc, "Lat")
trend_grace <- ncvar_get(grace_nc, "trend") #note that there are NAs for land datat
err_grace <- ncvar_get(grace_nc, "uncert")

grace_bdv <- data.frame(trend = as.numeric(trend_grace), std = as.numeric(err_grace), lon= lon, lat = lat)
grace_bdv <- grace_bdv[order(grace_bdv$lat),]
coordinates(grace_bdv) <- c("lon", "lat")

## Load BDV's ICE6G-VM5 in ewh mm/yr
gia_ewh <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA_BDV/GIA_ewh_ICE6G"), sep = ",", head = TRUE)
gia_ewh$lat <- -gia_ewh$lat
names(gia_ewh) <- c("lat", "lon", "trend")
gia_ewh$trend <- gia_ewh$trend*1000


## Load the Ocean polygon (-180, 180)
Ocean <- readOGR(dsn = paste0(dd, "WP1-BHM/maps/ne_110m_ocean"), layer = "ne_110m_ocean")

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

## Separate ocean and land obs
grace_coords <- grace_bdv@coords
grace_coords[,1] <- ifelse(grace_coords[,1]> 180, grace_coords[,1]-360, grace_coords[,1])
goid <- unlist(over(Ocean, SpatialPoints(coords= grace_coords, proj4string = CRS(proj4string(Ocean))), returnList = TRUE))
glid <- (1:nrow(grace_coords))[-goid]


## Set up prior
############################################################################################################
## Priors for mass
# -- for the ocean
mu_rg1 <- 1500/6371 # mean of rho
v_rg1 <- 1 # vague variace for rho
mu_sg1 <- 9 # mean of sigma
v_sg1 <- 18^2 # vaue variance for sigma
# -- for the land
mu_rg2 <- 600/6371 # mean of rho
v_rg2 <- 1 # vague variace for rho
mu_sg2 <- 16 # mean of sigma
v_sg2 <- 32^2 # vaue variance for sigma

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

## Generate mesh and build spde
############################################################################################################
## Genereate Fibonacci points on the sphere
fibo_points <- fiboSphere(N = 3e4, LL = FALSE)
mesh0 <- inla.mesh.2d(loc = fibo_points, cutoff = 0.01, max.edge = 1)
summary(mesh0)
mesh0ll <- Lxyz2ll(list(x = mesh0$loc[,1], y = mesh0$loc[,2], z = mesh0$loc[,3]))
mesh0ll <- cbind(mesh0ll$lon, mesh0ll$lat)
mesh0 <- dt.mesh.addon.posTri(mesh = mesh0, globe = TRUE) 
Tlonlat <- Lxyz2ll(list(x = mesh0$posTri[,1], y = mesh0$posTri[,2], z = mesh0$posTri[,3]))
mesh0$Trill <- cbind(lon = Tlonlat$lon, lat =Tlonlat$lat)
TinOcean <- unlist(over(Ocean, SpatialPoints(coords=mesh0$Trill, proj4string = CRS(proj4string(Ocean))), returnList=T))
TAll <- 1:mesh0$t
ToutOcean <- TAll[-TinOcean]
Omega = dt.Omega(list(TinOcean, 1:mesh0$t), mesh0)


mesh_ocean <- mesh.sub(mesh0, Omega, 1)
mesh_land <- mesh.sub(mesh0, Omega, 2)
summary(mesh_ocean)
summary(mesh_land)

## mass
massQ <- pp.create.Q(meshes = list(mesh_ocean, mesh_land), initial.theta = c(lrhog1, lrhog2, lsigmag1, lsigmag2))
## Transform the parameters for the SPDE_GMRF approximation
prior <- list(sigma = rbind(tsigmag1, tsigmag2), rho = rbind(trhog1, trhog2))
log.prior <- pp.create.prior.log.norm(prior.param = prior) 
mass_spde <- pp.inla.model(Q = massQ, log.prior=log.prior)


## Link process and data
############################################################################################################
## The data
ygrace_Ocean <- grace_bdv$trendgia[goid]
ygrace_Land <- grace_bdv$trendgia[-goid]

## Link GRACE to mass
grace_locO <- do.call(cbind, Lll2xyz(lat = grace_bdv@coords[goid,2], lon = grace_bdv@coords[goid,1]))
grace_locL <- do.call(cbind, Lll2xyz(lat = grace_bdv@coords[-goid,2], lon = grace_bdv@coords[-goid,1]))
grace_loc <- do.call(cbind, Lll2xyz(lat = grace_bdv@coords[,2], lon = grace_bdv@coords[,1]))
A_graceall <- inla.spde.make.A(mesh = mesh0, loc = grace_loc)
A_graceO <- inla.spde.make.A(mesh = mesh_ocean, loc = grace_locO)
A_graceL <- inla.spde.make.A(mesh = mesh_land, loc = grace_locL)
A_grace_m <- bdiag(A_graceO, A_graceL)

## The error
prec_scale <- c(1/grace_bdv$std[goid]^2, 1/grace_bdv$std[-goid]^2)


## Create the mass and steric stacks
stkmass <- inla.stack(data = list(y=c(ygrace_Ocean, ygrace_Land)), A = list(A_grace_m),
                      effects = list(mass = 1:(mesh_ocean$n + mesh_land$n)), 
                      remove.unused = FALSE, tag = "mass")

## run inla
############################################################################################################
## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))

## The formular -- we add the constraint that mass change sum to zero
formula = y ~ -1 + f(mass, model = mass_spde,
                     extraconstr = list(A = matrix(1, nrow = 1, ncol = (mesh_ocean$n + mesh_land$n)), e = 0))

res_inla <- inla(formula, data = inla.stack.data(stkmass), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stkmass), compute =TRUE), verbose = TRUE)