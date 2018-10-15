###################################################################
##                    GIA Data Screening                         ##
## This script read in and do screening of the GIA data produced ##
## from the ICE-6G forward model.                                ##
## We will use this data to:                                     ##
## 1. Generate a mesh with similar resolution                    ##
## 2. Use the data as prior mean of the GP(GMRF)                 ##
## 3. Learn initial values for lenght scale and variance         ##
###################################################################
library(rgl)
library(GEOmap)
library(INLA)
library(MVST)
library(fields)
library(gstat)
library(geoR)
#### 0 Read in GIA data
###################################################################
GIA_ice6g <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GIA_truth.txt", header = T)
polycoords <- GIA_ice6g[,c(6:13, 6,7)]
plist <- lapply(GIA_ice6g$ID, function(x) Polygons(list(Polygon(cbind(lon = as.numeric(polycoords[x, c(1,3,5,7,9)]), 
                                                                      lat = as.numeric(polycoords[x, c(2,4,6,8,10)])))), ID = x))
                                      
Plist <- SpatialPolygons(plist, proj4string = CRS("+proj=longlat"))
GIA_ice6g_sp <- SpatialPolygonsDataFrame(Plist, data = GIA_ice6g[,c("trend", "std")])
save(Plist, GIA_ice6g_sp, file = "GIA_sp_info.RData")
## Treat the gridded data as point data sampled at the centre of the grid
## converte the LonLat to xyz coordinates
GIA_loc6 <- do.call(cbind, Lll2xyz(lat = GIA_ice6g$y_center, lon = GIA_ice6g$x_center))
GIA_mu6 <- GIA_ice6g$trend

#### 1 Genreate a mesh by from the given points
###################################################################
## To get 1 by 1 degree resolution, we try to generate about 64800 triangles 
Mesh_GIA <- inla.mesh.2d(loc = GIA_loc6, cutoff = 0.1, max.edge = 0.38) ## small
Mesh_GIA <- inla.mesh.2d(loc = GIA_loc6, cutoff = 0.04, max.edge = 0.19) ## Medium
Mesh_GIA <- inla.mesh.2d(loc = GIA_loc6, cutoff = 0.0167, max.edge = 0.15) ## Large

GIA_spde <- inla.spde2.matern(Mesh_GIA)
summary(Mesh_GIA)
plot(Mesh_GIA, rgl = T)

MlocLL <- Lxyz2ll(list(x=Mesh_GIA$loc[,1], y = Mesh_GIA$loc[,2], z = Mesh_GIA$loc[,3]))
MlocLL$lon <- ifelse(MlocLL$lon <0, MlocLL$lon + 360, MlocLL$lon)
M_sp <- SpatialPoints(data.frame(lon = MlocLL$lon, lat = MlocLL$lat), proj4string = CRS("+proj=longlat")) #This convert GIA_ice6g a SpatialPointDataFrame
Midx <- over(M_sp, Plist)
GIA_muM <- GIA_mu6[Midx]

Mesh_GIA_sp <- SpatialPointsDataFrame(M_sp, data.frame(GIA_m = GIA_muM))


## Save all initial built up objects
save(GIA_ice6g_sp, Mesh_GIA, Mesh_GIA_sp, GIA_spde, file = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIAl.RData")


## Plot the GIA prior
open3d()
par3d(windowRect = c(100, 100, 900, 900))
layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
par3d(zoom = 0.8)
t_lim <- range(GIA_ice6g$trend)
t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
t_Cpal <- topo.colors(t_Clens, alpha=0)
t_Cols<- t_Cpal[round((GIA_ice6g$trend - t_lim[1])*100) + 1]
plot3d(GIA_loc6, col= t_Cols)

## plot the color bar
next3d(reuse = FALSE)
bgplot3d({z=matrix(1:t_Clens, nrow = t_Clens)
y=1
x=seq(t_lim[1],t_lim[2],len = t_Clens)
par(cex = 1.5, fin = c(8, 1), mai = c(0,0, 0.5, 0), oma = c(1, 0, 0, 0))
image(x,y,z,col = topo.colors(t_Clens),axes=FALSE,xlab="",ylab="")
title("The ICE6G GIA")
axis(1)})

writeWebGL(dir = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a", 
           filename= "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/GIA_globe.html",  # save the plot as an html
           width = 1000, reuse = TRUE)

#### 2 Set up GIA priors
###################################################################
## Use the GIA data to learn initial values for the covariance parameters
## Parameters needed for setting up the Matern covariance function
## nu -- mean square differentialbilty of the process -- poorly identified -- fix at 1 or 2
## var -- marginal variance  -- use the sample estimated sill and assuming no nugget
## kappa -- length scale -- range rho = sqrt(8nu/kappa) for correlation 0.1
## These can be determined by experts or data or both.

## plot and calculate the sample variogram from longlat coords
GIAdata <- data.frame(trend = GIA_ice6g$trend, lat = GIA_ice6g$y_center, long = GIA_ice6g$x_center)
coordinates(GIAdata) <- c("long", "lat")
proj4string(GIAdata) <- CRS("+proj=longlat")
## calculate the great circle distance from long lat so set projection to be FALSE
v0 <- variogram(trend ~ 1, data = GIAdata)
plot(v0)
## From the plot, sill = sigma = 5, range = 500. choose kappa = 1 nugget = 0 (assume)
f0 <- fit.variogram(v0, vgm(5, model = "Mat", kappa = 1, range = 500))
summary(f0)  ## fitted range = 343, convert to unit ball assumming earth radius 6371km = 0.0538
plot(v0, f0)

GIAgeo <- as.geodata(GIAdata)
vv0 <- likfit(GIAgeo, ini.cov.pars = c(5, 500), kappa = 1, fix.nugget = TRUE, nugget = 0, 
              cov.model = "matern", lik.method = "REML")

## Build the GMRF precision matrix for GIA
## if inla.mesh.fem which calls inla.fmesher.smorg detects S2, we can use
GIA_fem <- inla.mesh.fem(Mesh_GIA, order = 2)
Q_GIA <- Prec_from_SPDE_wrapper(M=GIA_fem$c1, K = GIA_fem$g1, nu = 1, desired_prec = 1/4.72, l = 0.0538)
## if not, we can only use inla function to build Q
## Setting pars in inla.matern
# d = 2 for S2
# alpha = nu + d/2, fractional operator, 0< alpha < 2, choose nu = 1, alpha = 1
# range -- initial set to be minimum of the mesh grid size
# log (tau) = theta1
# log (kappa) = theta2
#size0 <- min(c(diff(range(Mesh_GIA$loc[, 1])), diff(range(Mesh_GIA$loc[, 2])), diff(range(Mesh_GIA$loc[, 3]))))
range0 <- 0.0538
sigma0 <- sqrt(4.72)
kappa0 <- sqrt(8)/range0
tau0 <- 1/(4*pi*kappa0^2*sigma0^2)
GIA_Q0 <- inla.spde.precision(GIA_spde, theta=c(log(sqrt(tau0)), log(kappa0)))
GIA_Q1 <- inla.spde.precision(GIA_spde, theta=c(-1, 1.7))
xx <- inla.qsample(Q=GIA_Q1)
