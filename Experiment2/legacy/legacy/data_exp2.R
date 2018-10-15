#### Experiment 2c: Data processing overheads

## Uncomment when running this script alone
# if(Sys.info()["nodename"] == "IT034844"){
#   ## BS
#   wd <- "C:/ZSwork/"
#   dd <- "Z:/"
# }else if(Sys.info()["nodename"] == "IT064613"){
#   ## Maths
#   wd <-"/home/zs16444/coding/"
#   dd <-"/home/zs16444/globalmass/"
# }else{
#   ## server  
#   wd <- "~/"
#   dd <- "/./projects/GlobalMass/"
# }
# 
# library(rgdal); library(sp);library(GEOmap)

#########################################################################################################
#### Load all the data and make them into the desired format

### The altimetry data 
## -- Save alt_data, alt_data2 and alt_loc
alt_nc <- nc_open(paste0(dd,"WP3-Ocean/BHMinputs/SSH/William_Llovel/trend_SSH_CCI_200501_201512.nc"))
lon <- ncvar_get(alt_nc, "lon")
lat <- ncvar_get(alt_nc, "lat")
trend_ssh <- ncvar_get(alt_nc, "trend_ssh") #note that there are NAs for land datat
err_ssh <- ncvar_get(alt_nc, "err")
alt_data <- data.frame(trend_ssh = as.numeric(trend_ssh), err_ssh = as.numeric(err_ssh))
alt_data$lon <- rep(lon, 180) 
alt_data$lat <- rep(lat, each = 360)
alt_data2 <- na.omit(alt_data)
alt_loc <- do.call(cbind, Lll2xyz(lon = alt_data2$lon, lat = alt_data2$lat))

### Load the Land polygon (-180, 180)
## -- Save Land 
land0 <- readOGR(dsn = paste0(dd, "WP1-BHM/maps/ne_110m_land"), layer = "ne_110m_land")
## Remove small islands in the low lattidure(lat < 70) regrions (area < 300*300km)
areas <- sapply(land0@polygons, slot, "area")
lats <- sapply(land0@polygons, function(x){slot(x, "labpt")[2]})
smallareas <- areas <= 10
lowlats <- lats <= 70
lowsmall <- lowlats & smallareas
newids <- which(!lowsmall)
Land <- SpatialPolygons(land0@polygons[newids], proj4string = CRS(proj4string(land0)))


## The ICE6G-VM5 in vlm mm/yr
ice6g_vlm <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt"), header = T)
names(ice6g_vlm)[2:3] <- c("lon", "lat")
ice6g_vlm$lon <- ifelse(ice6g_vlm$lon >180, ice6g_vlm$lon - 360, ice6g_vlm$lon)
coordinates(ice6g_vlm) <-c("lon", "lat")
gridded(ice6g_vlm) <- TRUE
proj4string(ice6g_vlm) <- CRS(proj4string(Land))
ice6g_vlm$areas <- geosphere::areaPolygon(as(ice6g_vlm, "SpatialPolygons"))/(1000^2) # calculated the grid areas
Lid <- unlist(over(Land, as(ice6g_vlm, "SpatialPoints"), returnList = T)) # Find the grid in land
vlm <- sum(ice6g_vlm$trend[-Lid]*ice6g_vlm$areas[-Lid]/sum(ice6g_vlm$areas[-Lid])) # the average vlm in the ocean
print(paste("The vlm adjustment is ", vlm, "mm/yr"))

## Adjust the altimetry data by vlm
alt_data2$trend_vlm <- alt_data2$trend_ssh - vlm

### Load BDV's GRACE data 
## -- note the data are not evenly spaced!
## -- Save grace_sp, grace_df and grace_loc
grace_nc <- nc_open(paste0(dd, "WP2-SolidEarth/Bramha/BHMinput/GRACE_trends_Oceanadded_Bramha.nc"))
print(grace_nc)
lon <- ncvar_get(grace_nc, "Long")
lat <- ncvar_get(grace_nc, "Lat")
trend_grace <- ncvar_get(grace_nc, "trend") #note that there are NAs for land datat
err_grace <- ncvar_get(grace_nc, "uncert")

grace_df <- data.frame(trend = as.numeric(trend_grace), std = as.numeric(err_grace), lon= lon, lat = lat)
grace_sp <- grace_df
coordinates(grace_sp) <- c("lon", "lat")

### Load BDV's ICE6G-VM5 in ewh mm/yr
## -- Save gia_ewh
gia_ewh <- read.table(paste0(dd, "WP2-SolidEarth/BHMinputs/GIA_BDV/GIA_ewh_ICE6G"), sep = ",", head = TRUE)
gia_ewh$lat <- -gia_ewh$lat
names(gia_ewh) <- c("lat", "lon", "trend")
gia_ewh$trend <- gia_ewh$trend*1000
gia_sp <- gia_ewh
coordinates(gia_sp) <- c("lon", "lat")
gridded(gia_sp) <- TRUE
gia_sp$areas <- geosphere::areaPolygon(as(gia_sp, "SpatialPolygons"))/(1000^2)

## Remove GIA from GRACE -- use grace - gia grid value which grace falls in
proj4string(gia_sp)<- proj4string(grace_sp) <- CRS("+proj=longlat")
gridded(gia_sp) <- TRUE
giaG <- over(grace_sp, as(gia_sp, "SpatialPolygons"))
grace_df$trendgia <- grace_sp$trendgia <- grace_sp$trend - gia_sp$trend[giaG]
grace_loc <- do.call(cbind, Lll2xyz(lon = grace_df$lon, lat = grace_df$lat))

### Up to now, I want to save (alt_data, alt_data2, alt_loc, grace_sp, grace_df, grace_loc, Land)

#########################################################################################################
#### Some data check -- uncomment if want to check again
## GRACE sum to zero -- equal area -- simple average
# mean(grace_bdv$trend)
# 
# ## GIA_ewh sum to zero
# sum(gia_sp$trend*gia_sp$areas)/sum(gia_sp$areas)
# 
# ## GRACE adjusted by GIA also sum to zero
# mean(grace_bdv$trendgia)

#########################################################################################################

#### Create buffer zones and independent zones

### Create buffer zones for large land polygons.
## -- The land polygons are devided into large and small groups.
## -- Assuming the correlation lenght on the land can be as large as 600km, 
## -- then "Small" is defined by area < 400,000 km^2. This also cover most of small islands.
## -- The rest are "large", and we remove the polygon hole that represent the Caspian Sea.

### Buffer Description
## -- We create two buffer zones, one inside the land polygons and the other outside.
## -- Any points inside the buffer zone are assumed to be independent spatially.
## -- Buffer zone width: inner = 200km, outer = 400km

### 1 Remove the Caspian Sea
holes <- which(sapply(Land@polygons, function(x) length(x@Polygons) > 1))
caspian <- no_caspian <- Land@polygons[[holes]]
no_caspian@Polygons <- no_caspian@Polygons[1]
no_caspian@plotOrder <- 1L

## The land SpatialPolygons with no Caspian sea
Land2 <- Land
Land2@polygons[[holes]] <- no_caspian 
slot(Land2, "polygons") <- lapply(slot(Land2, "polygons"), "comment<-", NULL) 

## The Polygons list of Caspian Sea
caspian <- caspian@Polygons[2]

### 2 Identify the small polygons
smallid <- which(sapply(Land2@polygons, slot, "area") < 400)
Spoly <- Land2@polygons[smallid]
Lpoly <- Land2@polygons[-smallid]

### 3 Create the inner buffer and outer buffer for the large polygons
## All the polygons except Antarctica
p1 <- SpatialPolygons(Lpoly[-1], proj4string = CRS(proj4string(Land)))
## Transform the CRS of the large polygons for the bufferring operation
## Use a proper CRS: epsg:3857, see https://epsg.io/3857
## Works for lats within 85 -- need another one for Antarctica
p1 <- spTransform(p1, CRS("+init=epsg:3857"))
bp1in <- gBuffer(p1, width = -200*1000 )
bp1out <- gBuffer(p1, width = 400*1000 )
## Remove small area regions and holes in the buffers
smallid <- which(sapply(bp1in@polygons[[1]]@Polygons, slot, "area") < 400*1000*1e6)
bp1in <- SpatialPolygons(list(Polygons(bp1in@polygons[[1]]@Polygons[-smallid], ID = 1)),proj4string = CRS("+init=epsg:3857") )

holesid <- which(sapply(bp1out@polygons[[1]]@Polygons, slot, "hole"))
bp1out <- SpatialPolygons(list(Polygons(bp1out@polygons[[1]]@Polygons[-holesid], ID = 1)),proj4string = CRS("+init=epsg:3857") )

## Use another one for Antarctica, see https://epsg.io/5481
p2 <- SpatialPolygons(Lpoly[1], proj4string = CRS(proj4string(Land)))
p2 <- spTransform(p2, CRS("+init=epsg:5481"))
bp2in <- gBuffer(p2, width = -200*1000 )
bp2out <- gBuffer(p2, width = 400*1000 )
## !!!Note: no need to project them back! Project the points when using over!!!

### Now create the indpendent zones
## -- The independent zones are defined to include the followings
## -- Small area regions discarded at the beginning of creating buffer zones
## -- The difference band regions between inner and outer buffers
## -- The Caspian sea 
## -- Note that these can not be stored in one SpationPolygon Objects due to the different CRS!!!

## Check any small area regions above 80 degree in lats -- no so we can use the epsg:3857
highlats <- any(sapply(Spoly, function(x) abs(slot(x, "labpt")[2]) > 85))
## combine the caspian see and these regions into a srl
srl <- sapply(Spoly, function(x) x@Polygons[[1]])
sPoly <- Polygons(srl, ID = 1)
sPPoly <- SpatialPolygons(list(sPoly), proj4string = CRS(proj4string(Land)))
sPPoly <- spTransform(sPPoly, CRS("+init=epsg:3857"))

## Create the buffer band polygons (outer as polygon, inner as hole)
bufferbands1 <- gDifference(bp1out, bp1in, byid=TRUE) 

## Combine the small regions and the buffer bands 1
Indzone1 <- gUnion(bufferbands1, sPPoly)

## combine the caspian sea
caspian[[1]]@hole <- FALSE
caspian[[1]]@ringDir <- 1L
caspian <- SpatialPolygons(list(Polygons(caspian, ID = "caspian")), proj4string =CRS(proj4string(Land)))
caspian <- spTransform(caspian, CRS("+init=epsg:3857") )
Indzone1 <- gUnion(Indzone1, caspian)

## The second Indzone is the buffer bands2
Indzone2 <- bufferbands2 <- gDifference(bp2out, bp2in, byid=TRUE) 

plot(Indzone1, col = 2)
plot(Indzone2, col = 2)

## Create the Land and Ocean Polygons after removing the bufferzone
Land1 <- bp1in
Land2 <- bp2in

## We do not define the Ocean polygon explicilt here, since it is a tedious work
## We assume anything out side Indzones and Lands are Ocean

#########################################################################################################



#########################################################################################################
#### Priors setting

## Priors for mass
# -- Beta distribution for the range
# -- set the mode and let beta = 2*alpha
# -- mode = (alpha - 1)/(alpha + beta -2)
m_rg1 <- 1200/6371 # mode of rho_ocean
m_rg2 <- 600/6371  # mode of rho_land
# -- transform to beta parameters

Tbeta <- function(mode, beta.s){
  alpha <- (2*mode-1)/(mode*(1+beta.s)-1)
  beta <- beta.s*alpha
  return(c(alpha, beta))
}

trhog1 <- Tbeta(mode = m_rg1, beta.s = 3)
trhog2 <- Tbeta(mode = m_rg2, beta.s = 6)
## Plot the prior
xx <- seq(0,1, 0.01)
yy1 <- dbeta(xx, trhog1[1], trhog1[2])
yy2 <- dbeta(xx, trhog2[1], trhog2[2])
plot(xx, yy1, type = "l")
plot(xx, yy2, type = "l")
# -- initial values for rho_ocean and rho_land
lrhog1 <- log(m_rg1)
lrhog2 <- log(m_rg2)

## -- Gamma distribution for the variance
## -- set mean and variance
## -- mean = alpha/beta, variance = alpha/beta^2
m_sg1 <- 9        # mean of sigma_ocean
v_sg1 <- (9)^2/2 # variance of sigma_ocean
m_sg2 <- 19       # mean of sigma_land
v_sg2 <- (19)^2/2 # variance of sigma_land

TGamma <- function(mean, variance){
  beta <- mean/variance
  alpha <- mean^2/variance
  return(c(alpha, beta))
}

tsigmag1 <- TGamma(mean = m_sg1, variance = v_sg1)
tsigmag2 <- TGamma(mean = m_sg2, variance = v_sg2)
xx <- seq(0, 50, 1)
yy3 <- dgamma(xx, shape = tsigmag1[1], rate = tsigmag1[2])
yy4 <- dgamma(xx, shape = tsigmag2[1], rate = tsigmag2[2])
plot(xx, yy3, type = "l")
plot(xx, yy4, type = "l")

#-- initial values for sigma_Ocean and sigma_land
lsigmag1 <- log(m_sg1)
lsigmag2 <- log(m_sg2)

## Priors for steric
mu_r2 <- 2000/6371 # mean of rho
v_r2 <- 1 # vague variace for rho
mu_s2 <- 8 # mean of sigma
v_s2 <- 16^2 # vaue variance for sigma

## -- steric
trho2 <- Tlognorm(mu_r2, v_r2)
tsigma2 <- Tlognorm(mu_s2, v_s2)

lsigma2 <- tsigma2[1]
theta2_ss <- tsigma2[2]
lrho2 <- trho2[1]
theta2_rs <- trho2[2]
lkappa2 <- log(8)/2 - lrho2
ltau2 <- 0.5*log(1/(4*pi)) - lsigma2 - lkappa2


# ## create the ocean mesh
# ## Load the Land polygon (-180, 180)
# Land <- readOGR(dsn = paste0(dd, "WP1-BHM/maps/ne_110m_land"), layer = "ne_110m_land")
# ## Remove small islands in the low lattidure(lat < 70) regrions (area < 300*300km)
# areas <- sapply(Land@polygons, slot, "area")
# lats <- sapply(Land@polygons, function(x){slot(x, "labpt")[2]})
# smallareas <- areas <= 10
# lowlats <- lats <= 70
# lowsmall <- lowlats & smallareas
# newids <- which(!lowsmall)
# newLand <- SpatialPolygons(Land@polygons[newids], proj4string = CRS(proj4string(Land)))
# 
# TinLand <- unlist(over(newLand, SpatialPoints(coords=mesh0$Trill, proj4string = CRS(proj4string(newLand))), returnList=T))
# TAll <- 1:mesh0$t
# TinOcean <- TAll[-TinLand]
# Omega2 = dt.Omega(list(TinOcean, 1:mesh0$t), mesh0)
# 
# mesh_ocean <- mesh.sub(mesh0, Omega2, 1)

