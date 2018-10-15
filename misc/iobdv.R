## Load MS's GRACE data
gia_bdv <- read.table("c:/ZSwork/experiments/GIA_ewh_ICE6G", sep = ",", head = TRUE)
gia_bdv$lat <- -gia_bdv$lat
coordinates(gia_bdv) <- c("long", "lat")
gridded(gia_bdv) <- TRUE
gia_bdv$areas <- geosphere::areaPolygon(as(gia_bdv, "SpatialPolygons"))/(1000^2)
sum(gia_bdv$ewh.meters.*gia_bdv$areas)/sum(gia_bdv$areas)
spplot(gia_bdv, "ewh.meters.")

grace_loc <-  read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GRACE/GRACE_v02.3b_loc_v03.txt"), skip = 1)
n_grace <- ncol(grace_loc)

grace_nc <- nc_open("c:/ZSwork/experiments/GRACE_trends_Bramha.nc")
print(grace_nc)
lon <- ncvar_get(grace_nc, "Long")
lat <- ncvar_get(grace_nc, "Lat")
trend_grace <- ncvar_get(grace_nc, "trend") #note that there are NAs for land datat
err_grace <- ncvar_get(grace_nc, "uncert")

grace_bdv <- data.frame(trend = as.numeric(trend_grace), std = as.numeric(err_grace), lon= lon, lat = lat)
coordinates(grace_bdv) <- c("lon", "lat")

grace_bdvg <- grace_bdv
gridded(grace_bdvg) <- TRUE
spplot(grace_bdvg, "trend")
