## This script convert the following to ncdf fiels
## -- ice6g GIA data (vertical motion)
## -- ice6g GIA data (water equivalence height)
## -- GRACE super mascons data
## -- BHM predicted SSh
## -- BHM predicted Mass


## Math
# wd <-"/home/zs16444/"
# dd <-"/home/zs16444/globalmass/"

## BS
wd <- "C:/ZSwork/"
dd <- "Z:/" 

## load library and souce functions
library(raster)
library(ncdf4)
library(sp)
source(paste0(wd,"gmrcode/BHM_sphere/ncdfio.R"))


## ice6g vertical motion
ice6g <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt"), header = T)
ice6g$x_center <- ifelse(ice6g$x_center < 0, ice6g$x_center+360, ice6g$x_center)
ice6g2<- ice6g[order(ice6g$y_center,ice6g$x_center ), 2:5] # re-order the data according to the coordinates
names(ice6g2) <- c("lon", "lat", "trend", "std")

df <- ice6g2
fname <- "ice6g_VM5_vlm.nc"
vars <- list(c("trend", "trend", "vertical land motion trend", "mm/yr"), 
             c("std", "std", "standard error of vertical land motion", "mm/yr"))
title <- "ICE6g-VM5 GIA solution (vertical land motion)"

df2ncdf(df, fname, vars, title = title, append = FALSE)


## ice6g water equivalence
ice6g <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5_ewh.txt"), header = T)
ice6g$x_center <- ifelse(ice6g$x_center < 0, ice6g$x_center+360, ice6g$x_center)
ice6g2<- ice6g[order(ice6g$y_center,ice6g$x_center ), 2:5] # re-order the data according to the coordinates
names(ice6g2) <- c("lon", "lat", "trend", "std")

df <- ice6g2
fname <- "ice6g_VM5_ewh.nc"
vars <- list(c("trend", "trend", "water equivalence height trend", "mm/yr"), 
             c("std", "std", "standard error of water equivalence height motion", "mm/yr"))
title <- "ICE6g-VM5 GIA solution (water equivalence height)"

df2ncdf(df, fname, vars, title = title, append = FALSE)

## GRACE super mascons
grace_data <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_trends_v03.txt"), header = T)
grace_loc <-  read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_loc_v03.txt"), skip = 1)
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
SpPolygon <- SpatialPolygons(Polygons_list, proj4string = CRS("+proj=longlat +epllps=WGS84"))

grace_sp <- SpatialPolygonsDataFrame(SpPolygon,grace_data)
grace_grid <- SpatialPoints(coords = cbind(ice6g2$lon, ice6g2$lat), proj4string = CRS(proj4string(grace_sp)))
grace_gv<- over(grace_grid, grace_sp)$mmweq
grace_df <- data.frame(lon = ice6g2$lon, lat = ice6g2$lat, trend = grace_gv)

fname <- "graceSupM.nc"
vars <- list(c("trend", "trend", "Grace super Mascon average trend", "mm/yr"))
title <- "GRACE super mascon trend"

df2ncdf(grace_df, fname, vars, title = title, append = FALSE)

## BHM predictd SSH
load(paste0(dd, "WP1-BHM/Experiment2a/exp2a_ssh.RData"))
SSH_pred <- res_SSH$SSH_pred

fname <- "BHM-SSH.nc"
vars <- list(c("mean", "trend", "BHM predicted SSH trend", "mm/yr"), 
             c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr"))
title <- "BHM-predicted SSH"

df2ncdf(SSH_pred, fname, vars, title = title, append = FALSE)

## BHM predictd Mass
load(paste0(dd, "WP1-BHM/Experiment2a/exp2a_M2.RData"))
M_pred <- res_M2$M_pred

fname <- "BHM_mass.nc"
vars <- list(c("mean", "trend", "BHM predicted mass trend", "mm/yr"), 
             c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr"))
title <- "BHM-predicted mass"

df2ncdf(M_pred, fname, vars, title = title, append = FALSE)

## BHM predictd Mass 2
load(paste0(dd, "WP1-BHM/Experiment2a/exp2a_M2.RData"))
M_pred <- res_M2$M_pred
M_pred$diff <- grace_sp@data$pred_diff

fname <- "BHM_mass.nc"
vars <- list(c("mean", "trend", "BHM predicted mass trend", "mm/yr"), 
             c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr"))
title <- "BHM-predicted mass"

df2ncdf(M_pred, fname, vars, title = title, append = FALSE)

