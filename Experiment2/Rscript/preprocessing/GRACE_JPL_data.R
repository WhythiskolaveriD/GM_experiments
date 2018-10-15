#### Assemble GRACE GSFC data ####
## This file assemble the GRACE GSFC data and final product includes
## (1) 1 degree equi-area data, with GIA removed trend, inflated error and inflatd error + GIA error
## (2) Catchment region data, with GIA removed trend, inflated error and inflated error + GIA error

## Load the equi-areal grace data and polygons
## Load grace data
library(ncdf4)
library(sp)
grace_nc <- nc_open(paste0(dd_root, "WP2-SolidEarth/Bramha/Datsets/GRACE/JPL/GRACE_trends_jpl_BDV.nc"))
print(grace_nc)
lon <- ncvar_get(grace_nc, "Long")
lat <- ncvar_get(grace_nc, "Lat")
trend_grace <- ncvar_get(grace_nc, "trend") #note that there are NAs for land datat
err_grace <- ncvar_get(grace_nc, "uncert")

grace_raw <- data.frame(trend = as.numeric(trend_grace), std = as.numeric(err_grace), lon= lon, lat = lat)
## note that the current JPL dataset has -180 to 180 longitude which is different from other data
## transform this to be 0 - 360 so that the processing will remain the same as the others
grace_raw$lon <- ifelse(grace_raw$lon < 0, grace_raw$lon + 360, grace_raw$lon) 
coordinates(grace_raw) <- c("lon", "lat")

## Load BDV's SPH GIA mm/yr
gia_ewh <- read.table(paste0(dd_root, "WP3-Ocean/BHMinputs/GIA/ice6g_grace_mmyr_fromSphHarBdV.txt"), sep = "\t", skip = 1)
names(gia_ewh) <- c("ID", "lat", "lon", "trend")
gia_ewh$lon <- ifelse(gia_ewh$lon < 0, gia_ewh$lon + 360, gia_ewh$lon)
gia_ewh$error <- abs(0.2*gia_ewh$trend) # assume the errors are 20% of the gia trend

## Remove GIA from GRACE -- use grace - gia grid value which grace falls in
gia_sp <- gia_ewh
coordinates(gia_sp) <- c("lon", "lat")
gridded(gia_sp) <- TRUE
gia_sp$areas <- geosphere::areaPolygon(as(gia_sp, "SpatialPolygons"))/(1000^2)

proj4string(gia_sp)<- proj4string(grace_raw) <- CRS("+proj=longlat")
giaG <- over(grace_raw, as(gia_sp, "SpatialPolygons"))
grace_raw$trendgia <- grace_raw$trend - gia_sp$trend[giaG]


## Load polygons
load(paste0(dd_root, "WP1-BHM/maps/GRACE/landcatchment_usefull.rda"))
load(paste0(dd_root, "WP1-BHM/maps/GRACE/oceanMas_poly.rda"))

## Build the mapping matrix from points to polygons
library(Matrix)
ALand <- Matrix(0, nrow = nrow(SPdf_useful@data), ncol = nrow(grace_raw@data))
AAnt <- Matrix(0, nrow = nrow(Antarctica2@data), ncol = nrow(grace_raw@data))
AGreen <- Matrix(0, nrow = nrow(Greenland2@data), ncol = nrow(grace_raw@data))
AOcean <- Matrix(0, nrow = nrow(Ocean_spf@data), ncol = nrow(grace_raw@data))

## convert the mesh vertice location to longlat
coords0 <- coords1 <- grace_raw@coords
coords0[,1] <- ifelse(coords0[,1] < 0, coords0[,1] + 360, coords0[,1])
coords1[,1] <- ifelse(coords0[,1] > 180, coords0[,1] - 360, coords0[,1])
gracell0 <- SpatialPoints(coords = coords0)
gracell1 <- SpatialPoints(coords = coords1)
proj4string(gracell0) <- proj4string(gracell1) <- proj4string(SPdf_useful)

## Land
landid <- over(gracell1, SPdf_useful)$DBID
IDs <- SPdf_useful$DBID
for (i in 1:(length(IDs))){
  idtemp <- which(landid == IDs[i])
  ALand[i, idtemp] <- 1/length(idtemp)
}

## Ocean
Oceanid <- over(gracell0, Ocean_spf)$ID
IDs2 <- Ocean_spf$ID
for (i in 1:(length(IDs2))){
  idtemp <- which(Oceanid == IDs2[i])
  AOcean[i, idtemp] <- 1/length(idtemp)
}

## Greenland
gracell1 <- spTransform(gracell1, CRS(proj4string(Greenland2)))
Greenid <- over(gracell1, Greenland2)$SUBREGION1
IDs3 <- Greenland2$SUBREGION1
for (i in 1:(length(IDs3))){
  idtemp <- which(Greenid == IDs3[i])
  AGreen[i, idtemp] <- 1/length(idtemp)
}

## Antarctica
gracell1 <- spTransform(gracell1, CRS(proj4string(Antarctica2)))
Antid <- over(gracell1, Antarctica2)$Subregion
IDs4 <- Antarctica2$Subregion
for (i in 1:(length(IDs4))){
  idtemp <- which(Antid == IDs4[i])
  AAnt[i, idtemp] <- 1/length(idtemp)
}

Aall <- rbind(ALand, AGreen, AAnt, AOcean)


## Build the error covariance matrix under the assumption that r = 300km with matern kernel
## This is likely to be very slow because we need to calculate the distance between 40k+ points!
## So we use the SPDE approach to create an approximated precision matrix
## First build a mesh on globe 
source(paste0(wd_root,"gmrcode/BHM_sphere/functions.R"))
library(INLA)
fibo_points <- fiboSphere(N = 3e4, LL = FALSE)
mesh0 <- inla.mesh.2d(loc = fibo_points, cutoff = 0.01, max.edge = 1)

## Then create the precision matrix for mass
range0 <- 300/6371
sigma0 <- 1
lkappa0 <- log(8)/2 - log(range0)
ltau0 <- 0.5*log(1/(4*pi)) - log(sigma0) - lkappa0

spde <- inla.spde2.matern(mesh0, B.tau = matrix(c(ltau0, -1, 1),1,3),
                          B.kappa = matrix(c(lkappa0, 0, -1), 1,3), 
                          theta.prior.mean = c(0,0), theta.prior.prec = c(0.1, 1))
Q <- inla.spde.precision(spde, theta = c(0,0))
Matern_cov <- inla.qinv(Q)

## Then create the precision matrix for gia
range0 <- 200/6371
sigma0 <- 1
lkappa0 <- log(8)/2 - log(range0)
ltau0 <- 0.5*log(1/(4*pi)) - log(sigma0) - lkappa0

spde <- inla.spde2.matern(mesh0, B.tau = matrix(c(ltau0, -1, 1),1,3),
                          B.kappa = matrix(c(lkappa0, 0, -1), 1,3), 
                          theta.prior.mean = c(0,0), theta.prior.prec = c(0.1, 1))
Q <- inla.spde.precision(spde, theta = c(0,0))
Matern_cov2 <- inla.qinv(Q)


## The covariance for the GRACE raw data
## Convert coordinates of GRACE raw data to xyz
raw.xyz <- do.call(cbind, GEOmap::Lll2xyz(lat = grace_raw@coords[,2], lon = grace_raw@coords[,1]))
Araw <- inla.spde.make.A(mesh = mesh0, loc = raw.xyz)
raw_cor <- cov2cor(Araw %*% Matern_cov %*% t(Araw))

sd_raw <- grace_raw$std
sd_raw <- Diagonal(n = length(sd_raw), sd_raw)
raw_cov1 <- sd_raw %*% raw_cor %*% sd_raw

## The inflation for raw data
a.raw1 <- sum(raw_cov1)/sum(diag(raw_cov1))
raw_error1<- sqrt(diag(raw_cov1)*a.raw1)


## The covariance for the gia raw data
## Convert coordinates of gia raw data to xyz
raw.xyz <- do.call(cbind, GEOmap::Lll2xyz(lat = gia_sp@coords[giaG,2], lon = gia_sp@coords[giaG,1]))
Araw <- inla.spde.make.A(mesh = mesh0, loc = raw.xyz)
raw_cor <- cov2cor(Araw %*% Matern_cov2 %*% t(Araw))

sd_raw <- gia_sp$error[giaG]
sd_raw <- Diagonal(n = length(sd_raw), sd_raw)
raw_cov2 <- sd_raw %*% raw_cor %*% sd_raw

## The inflation for raw data with gia error
a.raw2 <- sum(raw_cov2)/sum(diag(raw_cov2))
raw_error2<- sqrt(diag(raw_cov2)*a.raw2 + raw_error1^2)

## The land correlation
land_cov1 <- ALand %*% raw_cov1 %*% t(ALand)
land_cov2 <- ALand %*% raw_cov2 %*% t(ALand)
## The Ocean correlation
ocean_cov1 <- AOcean %*% raw_cov1 %*% t(AOcean)
ocean_cov2 <- AOcean %*% raw_cov2 %*% t(AOcean)
## Greenland 
Green_cov1 <- AGreen %*% raw_cov1 %*% t(AGreen)
Green_cov2 <- AGreen %*% raw_cov2 %*% t(AGreen)
## Antarctica
Ant_cov1 <- AAnt %*% raw_cov1 %*% t(AAnt)
Ant_cov2 <- AAnt %*% raw_cov2 %*% t(AAnt)

## Now find the inflation for land
a.land1 <- sum(land_cov1)/sum(diag(land_cov1))
a.land2 <- sum(land_cov2)/sum(diag(land_cov2))
land_error1 <- sqrt(diag(land_cov1)*a.land1)
land_error2 <- sqrt(land_error1^2 + diag(land_cov2)*a.land2)
land_val <- drop(ALand %*% grace_raw$trendgia)

## For Greenland
a.green1 <- sum(Green_cov1)/sum(diag(Green_cov1))
a.green2 <- sum(Green_cov2)/sum(diag(Green_cov2))
green_error1 <- sqrt(diag(Green_cov1)*a.green1)
green_error2 <- sqrt(green_error1^2 + diag(Green_cov2)*a.green2)
green_val <- drop(AGreen %*% grace_raw$trendgia)

## For Antarctica
a.Ant1 <- sum(Ant_cov1)/sum(diag(Ant_cov1))
a.Ant2 <- sum(Ant_cov2)/sum(diag(Ant_cov2))
Ant_error1 <- sqrt(diag(Ant_cov1)*a.Ant1)
Ant_error2 <- sqrt(Ant_error1^2 + diag(Ant_cov2)*a.Ant2)
Ant_val <- drop(AAnt %*% grace_raw$trendgia)

## Now find the inflation for Ocean
a.ocean1 <- sum(ocean_cov1)/sum(diag(ocean_cov1))
a.ocean2 <- sum(ocean_cov2)/sum(diag(ocean_cov2))
ocean_error1 <- sqrt(diag(ocean_cov1)*a.ocean1)
ocean_error2 <- sqrt(ocean_error1^2 + diag(ocean_cov2)*a.ocean2)
ocean_val <- drop(AOcean %*% grace_raw$trendgia)

## Combine the two polygon dataframe to one
library(maptools)
SPdf_useful@data <- data.frame(area = SPdf_useful@data$area, trend = land_val, std = land_error1, std2 = land_error2)
Ocean_spf@data <- data.frame(area = Ocean_spf@data$area, trend=ocean_val, std=ocean_error1, std2 = ocean_error2)
Greenland<-Greenland2
Greenland@data <- data.frame(trend = green_val, std = green_error1, std2 = green_error2)
Antarctica <- Antarctica2
Antarctica@data <- data.frame(trend = Ant_val, std = Ant_error1, std2 = Ant_error2)

SPdf_useful <- spChFIDs(SPdf_useful, paste0("land", row.names(SPdf_useful), sep="."))
Ocean_spf <- spChFIDs(Ocean_spf, paste0("ocean", row.names(Ocean_spf), sep="."))
Land_spf <- SPdf_useful

grace_sp <- grace_raw
grace_sp$std1 <- raw_error1
grace_sp$std2 <- raw_error2

## Load the Ocean polygon (-180, 180)
Ocean <- readOGR(dsn = paste0(dd_root, "WP1-BHM/maps/ne_110m_ocean"), layer = "ne_110m_ocean")
Land <- readOGR(dsn = paste0(dd_root, "WP1-BHM/maps/ne_110m_land"), layer = "ne_110m_land")


save(Land_spf, Ocean_spf, Antarctica, Greenland, grace_sp, gia_sp, Ocean, Land,
     file = paste0(dd_root, "WP1-BHM/Experiment2/data/Exp2_JPL_grace.rda"))










# ################################################################################################################
# ##### Plot for sanity check
# ## Scale the data for plot
# SPdf_useful$trend2 <- ifelse(abs(SPdf_useful$trend) > 19, sign(SPdf_useful$trend)*20, SPdf_useful$trend)
# colpal1 <- rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac'))
# spplot(SPdf_useful, "trend2", at =c(-20, -15, -10, -5, 0, 5, 10, 15, 20), col.regions = colpal1)
# 
# SPdf_useful$stds <- ifelse(abs(SPdf_useful$std) > 4, 4.5, SPdf_useful$std)
# SPdf_useful$stds2 <- ifelse(abs(SPdf_useful$std2) > 4, 4.5, SPdf_useful$std)
# colpal2 <- c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#8c2d04')
# spplot(SPdf_useful, c("stds", "stds2"), at = c(0, 0.5, 1, 2, 3, 4, 5), col.regions = colpal2 )
# 
# spplot(Greenland, "trend", at =c(-161, -100, -75, -50, -30, -20, -10, -5, 0), col.regions = colpal1 )
# spplot(Greenland, "std", at = c(0, 0.5, 1, 1.5, 2, 2.5, 3.5), col.regions = colpal2 )
# spplot(Greenland, c("std", "std2"))
# 
# spplot(Antarctica, "trend", at =c(-200, -100, -50, -25, -10, 0, 10, 20, 35), col.regions = colpal1 )
# spplot(Antarctica, "std", at = c(0, 0.5, 1, 1.5, 2, 2.5, 3.5), col.regions = colpal2 )
# spplot(Antarctica, c("std", "std2"))
# 
# Ocean_spf$trend2 <- ifelse(abs(Ocean_spf$trend) > 19, sign(Ocean_spf$trend)*20, Ocean_spf$trend)
# spplot(Ocean_spf, "trend2", at =c(-20, -15, -10, -5, 0, 5, 10, 15, 20), col.regions = colpal1)
# 
# spplot(Ocean_spf, c("std2","std"))
# 
# 
