## Load data and priors
#########################################################################################################
## Try load from cache
cache_grace <- paste0(dd_root,"WP1-BHM/Experiment2/data/Exp2_GSFC_grace.rda")
cache_alt <- paste0(dd_root,"WP1-BHM/Experiment2/data/Exp2_GSFC_alt.rda")

# Note that these data are exactly the same as the GSFCgrid,
# since the data processing script save both grid and catchment region data
# This script is redundant but kept to be differenciate with the GSFCgrid script.
# This script can be deleted/merged in future.

## GRACE
if(file.exists(cache_grace)){
  load(cache_grace)
}else{
  source(paste0(wd_root, "experiments/Experiment2/Rscript/GRACE_GSFC_data.R")) 
  ## This script process and save both grid and catchment region data for GRACE GSFC data
}

## Altimetry
if(file.exists(cache_alt)){
  load(cache_alt)
}else{
  ## Load and process the data; then save as cache file
  ## The altimetry data
  alt_nc <- nc_open(paste0(dd_root,"WP3-Ocean/BHMinputs/SSH/trend_SSH_CCI_spatialGIA_OBD_2005_2015.nc"))
  lon <- ncvar_get(alt_nc, "lon")
  lat <- ncvar_get(alt_nc, "lat")
  trend_ssh <- ncvar_get(alt_nc, "trend_ssh") #note that there are NAs for land datat
  err_ssh <- ncvar_get(alt_nc, "err")
  nc_close(alt_nc)
  
  alt_data <- data.frame(trend_ssh = as.numeric(trend_ssh), err_ssh = as.numeric(err_ssh))
  alt_data$lon <- rep(lon, 180) 
  alt_data$lat <- rep(lat, each = 360)
  alt_data2 <- na.omit(alt_data)
  alt_loc <- do.call(cbind, Lll2xyz(lon = alt_data2$lon, lat = alt_data2$lat))
  
  save(alt_data, alt_data2, alt_loc, file = cache_alt)
}