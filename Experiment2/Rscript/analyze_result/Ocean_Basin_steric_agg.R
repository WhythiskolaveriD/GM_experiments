### Ocean basin aggregation
######################################################################################
### This script aggragat the steric trend at ocean basin level
library(ncdf4)
library(sp)


#### Load the baisn polygon and the steric data
######################################################################################
## The basins
load("Z:/WP1-BHM/maps/Ocean/basins.rda")

## function for loading the steric netcdf
nc2df <- function(file){
  ncfile <- nc_open(file)
  trend <- as.numeric(ncvar_get(ncfile, "trend"))
  error <- as.numeric(ncvar_get(ncfile, "trend_error"))
  lon <- ncvar_get(ncfile, "lon")
  lat <- ncvar_get(ncfile, "lat")
  nc_close(ncfile)
  nlon <- length(lon)
  nlat <- length(lat)
  df <- data.frame(lon = rep(lon, nlat), lat =  rep(lat, each = nlon),
                         trend = trend, error = error)
  return(df)
}

## Ensemble
ens_data <- nc2df(file = "Z:/WP3-Ocean/BHMinputs/steric/EnsembleOI_0_2000m_2005_2015.nc")
lattice::levelplot(error ~ lon + lat, data = ens_data)

## EN4
en4_data <- nc2df(file = "Z:/WP3-Ocean/BHMinputs/steric/EN4_4.2.1_0_2000m_2005_2015.nc")
lattice::levelplot(error ~ lon + lat, data = en4_data)

## IFREMER
IFR_data <- nc2df(file = "Z:/WP3-Ocean/BHMinputs/steric/IFREMER_ISAS15_0_2000m_2005_2015.nc")
lattice::levelplot(error ~ lon + lat, data = IFR_data)

## JAMSTEC
JAM_data <- nc2df("Z:/WP3-Ocean/BHMinputs/steric/JAMSTEC_MOAA_GPV_0_2000m_2005_2015.nc")
lattice::levelplot(error ~ lon + lat, data = JAM_data)

## Scripps
scrp_data <- nc2df("Z:/WP3-Ocean/BHMinputs/steric/Scripps_RG_ArgoClim_OI_0_2000m_2005_2015.nc")
lattice::levelplot(error ~ lon + lat, data = scrp_data)
scrp_data$trend[which(scrp_data$trend == 0)] <- NA ## The NA's here are all filled with value zero, somehow.
scrp_data$error[which(scrp_data$error == 0)] <- NA
lattice::levelplot(error ~ lon + lat, data = scrp_data)

### Aggregate the basin level trend and error
######################################################################################
## function for aggregation
spdata <- ens_data

sp_agg <- function(spdata, mod){
  ## data -- dataframe with lonlat defind point values
  ## mod -- the model name
  Basins <- as(allbasinsdf, "SpatialPolygons")
  ## Check the longitude of the data
  spdata$lon <- ifelse(spdata$lon < 0, spdata$lon + 360, spdata$lon)
  
  ## change the data to a SpatialPointDataFrame
  coordinates(spdata) <- c("lon", "lat")
  proj4string(spdata) <- CRS(proj4string(Basins)) 
  ids <- over(spdata, Basins)
  means <- sapply(1:6, function(x) mean(na.omit(spdata$trend[which(ids == x)])))
  sds <- sapply(1:6, function(x) sqrt(sum(na.omit(spdata$error[which(ids == x)])^2))/length(na.omit(spdata$error[which(ids == x)])))
  basindf <- data.frame(names = names(Basins), mean = means, sd = sds, model = mod)
  return(basindf)
}

ens_basin <- sp_agg(spdata = ens_data, mod = "ensemble")
en4_basin <- sp_agg(spdata = en4_data, mod = "EN4")
IFR_basin <- sp_agg(spdata = IFR_data, mod = "IFR")
JAM_basin <- sp_agg(spdata = JAM_data, mod = "JAM")
scrp_basin <- sp_agg(spdata = scrp_data, mod = "scripps")

basin_all <- do.call(rbind, list(ens_basin, en4_basin, IFR_basin, JAM_basin, scrp_basin))
write.csv(basin_all, file = "basin_mean.csv")

## convert the lattitude degree to radiant
deg2rad <- function(deg){
  rad <- abs(deg/180*pi)
  rad
}
## the area can be approximated by cos(radian(lat))


## Areal weighted average 
sp_agg <- function(spdata, mod, area = FALSE){
  ## data -- dataframe with lonlat defind point values
  ## mod -- the model name
  Basins <- as(allbasinsdf, "SpatialPolygons")
  ## Check the longitude of the data
  spdata$lon <- ifelse(spdata$lon < 0, spdata$lon + 360, spdata$lon)
  
  ## change the data to a SpatialPointDataFrame
  coordinates(spdata) <- c("lon", "lat")
  proj4string(spdata) <- CRS(proj4string(Basins)) 
  ids <- over(spdata, Basins)
  if(area){
    areas <- cos(deg2rad(spdata$lat))
    means <- sapply(1:6, function(x) sum(na.omit(spdata$trend[which(ids == x)]*
                                                   areas[which(ids == x)]/sum(areas[which(ids == x)]))))
    mean_global <- sum(na.omit(spdata$trend*areas/sum(areas)))
    sds <- sapply(1:6, 
                  function(x) sqrt(sum(na.omit(spdata$error[which(ids == x)]^2*
                                         areas[which(ids == x)]/(sum(areas[which(ids == x)]))^2))))
    sd_global <- sqrt(sum(na.omit(spdata$error^2*(areas/sum(areas))^2)))
  }
  else{
    means <- sapply(1:6, function(x) mean(na.omit(spdata$trend[which(ids == x)])))
    mean_global <- mean(na.omit(spdata$trend))
    sds <- sapply(1:6, 
                function(x) sqrt(sum(na.omit(spdata$error[which(ids == x)])^2))/
                  length(na.omit(spdata$error[which(ids == x)])))
    sd_global <- sqrt(sum(na.omit(spdata$error^2)))/length(na.omit(spdata$error))
    
  }
  basindf <- data.frame(names = c(names(Basins), "global"), mean = c(means, mean_global), 
                        sd = c(sds, sd_global), model = mod)
  return(basindf)
}

ens_basin <- sp_agg(spdata = ens_data, mod = "ensemble", area = TRUE)
en4_basin <- sp_agg(spdata = en4_data, mod = "EN4", area = TRUE)
IFR_basin <- sp_agg(spdata = IFR_data, mod = "IFR", area = TRUE)
JAM_basin <- sp_agg(spdata = JAM_data, mod = "JAM",area = TRUE)
scrp_basin <- sp_agg(spdata = scrp_data, mod = "scripps", area = TRUE)

basin_all2 <- do.call(rbind, list(ens_basin, en4_basin, IFR_basin, JAM_basin, scrp_basin))
write.csv(basin_all2, file = "basin_mean_weighted.csv")
