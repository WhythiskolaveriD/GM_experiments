## Load data and priors
#########################################################################################################
## Try load from cache
cache_data <- paste0(dd,"WP1-BHM/Experiment2c/data/exp2c_data.rda")
cache_vgram <- paste0(dd,"WP1-BHM/Experiment2c/data/exp2c_vgram.rda")
cache_mesh <- paste0(dd,"WP1-BHM/Experiment2c/data/exp2c_mesh.rda")

if(file.exists(cache_data)){
  load(cache_data)
}else{
  ## Load and process the data; then save as cache file
  ## The altimetry data
  alt_nc <- nc_open(paste0(dd,"WP3-Ocean/BHMinputs/SSH/trend_SSH_CCI_200501_201512.nc"))
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
  
  ## Load BDV's GRACE data -- note the data are not evenly spaced!
  grace_nc <- nc_open(paste0(dd, "WP2-SolidEarth/Bramha/BHMinput/GRACE_trends_Oceanadded_Bramha.nc"))
  print(grace_nc)
  lon <- ncvar_get(grace_nc, "Long")
  lat <- ncvar_get(grace_nc, "Lat")
  trend_grace <- ncvar_get(grace_nc, "trend") #note that there are NAs for land datat
  err_grace <- ncvar_get(grace_nc, "uncert")
  nc_close(grace_nc)
  
  grace_sp <- data.frame(trend = as.numeric(trend_grace), std = as.numeric(err_grace), lon= lon, lat = lat)
  coordinates(grace_sp) <- c("lon", "lat")
  
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
  
  ## Adjust the altimetry data by vlm
  alt_data2$trend_vlm <- alt_data2$trend_ssh - vlm
  
  ## GIA_ewh sum to zero
  gia_sp <- gia_ewh
  coordinates(gia_sp) <- c("lon", "lat")
  gridded(gia_sp) <- TRUE
  gia_sp$areas <- geosphere::areaPolygon(as(gia_sp, "SpatialPolygons"))/(1000^2)
  
  ## Remove GIA from GRACE -- use grace - gia grid value which grace falls in
  proj4string(gia_sp)<- proj4string(grace_sp) <- CRS("+proj=longlat")
  gridded(gia_sp) <- TRUE
  giaG <- over(grace_sp, as(gia_sp, "SpatialPolygons"))
  grace_sp$trendgia <- grace_sp$trend - gia_sp$trend[giaG]
  
  ## variagrams for setting up priors
  grace_coords <- grace_sp@coords
  grace_coords[,1] <- ifelse(grace_coords[,1]> 180, grace_coords[,1]-360, grace_coords[,1])
  goid <- unlist(over(Ocean, SpatialPoints(coords= grace_coords, proj4string = CRS(proj4string(Ocean))), returnList = TRUE))
  glid <- (1:nrow(grace_coords))[-goid]
  v_ocean <- gstat::variogram(trendgia ~ 1, data = grace_sp[goid,])
  v_land <- gstat::variogram(trendgia ~ 1, data = grace_sp[-goid,])
  
  save(alt_data, alt_data2, alt_loc, grace_sp, gia_sp, ice6g_vlm, vlm, Ocean,  file = cache_data)
}

## The variogram
#########################################################################################################
if(file.exists(cache_vgram)){
  load(cache_vgram)
}else{
  ## variagrams for setting up priors
  grace_coords <- grace_sp@coords
  grace_coords[,1] <- ifelse(grace_coords[,1]> 180, grace_coords[,1]-360, grace_coords[,1])
  goid <- unlist(over(Ocean, SpatialPoints(coords= grace_coords, proj4string = CRS(proj4string(Ocean))), returnList = TRUE))
  glid <- (1:nrow(grace_coords))[-goid]
  v_ocean <- gstat::variogram(trendgia ~ 1, data = grace_sp[goid,])
  v_land <- gstat::variogram(trendgia ~ 1, data = grace_sp[-goid,])
  save(goid, glid, v_ocean, v_land, file = cache_vgram)
}

## generate the mesh
#########################################################################################################
if(file.exists(cache_mesh)){
  load(cache_mesh)
}else{
  ## Genereate Fibonacci points on the sphere
  fibo_points <- fiboSphere(N = 3e4, LL = FALSE)
  mesh0 <- inla.mesh.2d(loc = fibo_points, cutoff = 0.01, max.edge = 1)
  summary(mesh0)
  
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

  save(mesh0, Omega, mesh_land, mesh_ocean, file = cache_mesh)
}


#### Set up the priors
#########################################################################################################
## Priors for mass
# -- for the ocean
mu_rg1 <- 1500/6371 # mean of rho
v_rg1 <- 1 # vague variace for rho
mu_sg1 <- 9 # mean of sigma
v_sg1 <- 18^2 # vaue variance for sigma
# -- for the land
mu_rg2 <- 700/6371 # mean of rho
v_rg2 <- 1 # vague variace for rho
mu_sg2 <- 16 # mean of sigma
v_sg2 <- 32^2 # vaue variance for sigma

## Priors for steric
mu_r2 <- 2000/6371 # mean of rho
v_r2 <- 1 # vague variace for rho
mu_s2 <- 15 # mean of sigma
v_s2 <- 30^2 # vaue variance for sigma

## Transform the parameters for a lognormal distribution
#########################################################################################################
## -- mass
trhog1 <- Tlognorm(mu_rg1, v_rg1)
tsigmag1 <- Tlognorm(mu_sg1, v_sg1)
trhog2 <- Tlognorm(mu_rg2, v_rg2)
tsigmag2 <- Tlognorm(mu_sg2, v_sg2)

lsigmag1 <- tsigmag1[1]
thetag1_ss <- tsigmag1[2]
lrhog1 <- trhog1[1]
thetag1_rs <- trhog1[2]
lkappag1 <- log(8)/2 - lrhog1
ltaug1 <- 0.5*log(1/(4*pi)) - lsigmag1 - lkappag1

lsigmag2 <- tsigmag2[1]
thetag2_ss <- tsigmag2[2]
lrhog2 <- trhog2[1]
thetag2_rs <- trhog2[2]
lkappag2 <- log(8)/2 - lrhog2
ltaug2 <- 0.5*log(1/(4*pi)) - lsigmag2 - lkappag2

## -- steric
trho2 <- Tlognorm(mu_r2, v_r2)
tsigma2 <- Tlognorm(mu_s2, v_s2)

lsigma2 <- tsigma2[1]
theta2_ss <- tsigma2[2]
lrho2 <- trho2[1]
theta2_rs <- trho2[2]
lkappa2 <- log(8)/2 - lrho2
ltau2 <- 0.5*log(1/(4*pi)) - lsigma2 - lkappa2