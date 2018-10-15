################################################################################
##### This script shows the effects of mesh on building the spde precision #####
################################################################################
## We use 3 different sets of locations to generate the mesh with small,      ##
## medium and large size.                                                     ##
## For each mesh generate the precicion matrix with the same parameters,      ##
## Map the marginal stardard errors of each mesh support                      ##
################################################################################

## 0 if run the script independently ##
# wkdir <- getwd()
## Priors mean and variance for rho and sigma
# mu_r <- 500/6371
# v_r <- (1000/6371)^2
# mu_s <- 20
# v_s <- 40^2


#### 1 load libraries, GPS data and  GIA prior
#library(rgl)
library(INLA)
if(runserver){
  INLA:::inla.dynload.workaround()
}#Use this on server with old gcc library
library(fields)
## GPS
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
GPSX <- ifelse(GPS_obs$x_center > 180, GPS_obs$x_center-360, GPS_obs$x_center)
GPSY <- GPS_obs$y_center
GPS_loc <- inla.mesh.map(cbind(GPSX, GPSY), projection = "longlat")
## GIA
GIA_ice6g <- read.table("experimentBHM/GIA_truth.txt", header = T)
GIA_ice6g$x_center <- ifelse(GIA_ice6g$x_center > 180, GIA_ice6g$x_center-360, GIA_ice6g$x_center)
GIA_loc6 <- inla.mesh.map(cbind(GIA_ice6g$x_center, GIA_ice6g$y_center), projection = "longlat")
GIA_mu6 <- GIA_ice6g$trend
GIA_sd6 <- GIA_ice6g$std


#### 2 Setup for spde parameters
source("glbm/Experiment1a/Rscript/MVSTplus.R")
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]

lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0



#### 3.1 Genreate a mesh by from the given GIA grid
## Find index near equator
## Find index near poles
ind_eq <- which(GIA_ice6g$y_center == 0.5)
ind_po <- which(GIA_ice6g$y_center == 89.5)
GIA_dist1 <- dist(GIA_loc6[ind_eq,])
GIA_dist2 <- dist(GIA_loc6[ind_po,])
summary(GIA_dist1)
summary(GIA_dist2)
u_dist <- summary(GIA_dist1)[1]
Mesh_GIAs <- list()
Mesh_GIAs[[1]] <- inla.mesh.2d(loc = GIA_loc6, cutoff = u_dist*8, max.edge = u_dist*12) ## small
Mesh_GIAs[[2]] <- inla.mesh.2d(loc = GIA_loc6, cutoff = u_dist*3, max.edge = u_dist*6) ## Medium
Mesh_GIAs[[3]] <- inla.mesh.2d(loc = GIA_loc6, cutoff = u_dist, max.edge = u_dist*3) ## Large
## Cutoff slightly larger than the grid width -- that's why some part of the generated grid is not regular
save(Mesh_GIAs, file = paste0(wkdir, "/experimentBHM/mesh_GIA.RData"))

pdf(file = paste0(wkdir, "/experimentBHM/prior_GIAmesh.pdf"), width = 8, height = 11)
mains <- paste("Prior marginal standard error -- ",  c("small mesh", "medium mesh", "large mesh"))
par(mfrow = c(3,1))
for (i in 1:3){
GIA_spde <- inla.spde2.matern(Mesh_GIAs[[i]], B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

Q00 <- inla.spde2.precision(GIA_spde, theta = c(0,0))
GIA_sprior <- sqrt(1/diag(Q00))
proj<- inla.mesh.projector(Mesh_GIAs[[i]], projection = "longlat", dims = c(361,181))
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_sprior)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = mains[i])
points(GPSX, GPSY, cex = 0.5)
}
dev.off()


#### 3.2 Genreate a mesh from a globe with equidistant points from equidistant lattice
mesh_regs <- list()
mesh_regs[[1]] <- inla.mesh.create(globe = 10)
mesh_regs[[2]] <- inla.mesh.create(globe = 20)
mesh_regs[[3]] <- inla.mesh.create(globe = 50)

save(mesh_regs, file = paste0(wkdir, "/experimentBHM/mesh_reg.RData"))

pdf(file = paste0(wkdir, "/experimentBHM/prior_Globemesh.pdf"), width = 8, height = 11)
mains <- paste("Prior marginal standard error -- ",  c("small mesh", "medium mesh", "large mesh"))
par(mfrow = c(3,1))
for (i in 1:3){
  GIA_spde <- inla.spde2.matern(mesh_regs[[i]], B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                                theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

  Q00 <- inla.spde2.precision(GIA_spde, theta = c(0,0))
  GIA_sprior <- sqrt(1/diag(Q00))
  proj<- inla.mesh.projector(mesh_regs[[i]], projection = "longlat", dims = c(361,181))
  image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_sprior)), col = topo.colors(40),
             xlab = "Longitude", ylab = "Latitude", main = mains[i])
  points(GPSX, GPSY, cex = 0.5)
}
dev.off()


#### 3.3 Genreate a mesh from GPS location
mesh_GPS <- list()
GPS_dist <- dist(GPS_loc)
summary(GPS_dist)
mesh_GPS[[1]] <- inla.mesh.2d(loc = GPS_loc, cutoff = u_dist*5,  max.edge = u_dist*12)
mesh_GPS[[2]] <- inla.mesh.2d(loc = GPS_loc, cutoff = u_dist*3, max.edge = u_dist*5)
mesh_GPS[[3]] <- inla.mesh.2d(loc = GPS_loc, cutoff = u_dist/3, max.edge = u_dist*2)
save(mesh_GPS, file = paste0(wkdir, "/experimentBHM/mesh_GPS.RData"))

pdf(file = paste0(wkdir, "/experimentBHM/prior_GPSmesh.pdf"), width = 8, height = 11)
mains <- paste("Prior marginal standard error -- ",  c("small mesh", "medium mesh", "large mesh"))
par(mfrow = c(3,1))
for (i in 1:3){
  GIA_spde <- inla.spde2.matern(mesh_GPS[[i]], B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                                theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))

  Q00 <- inla.spde2.precision(GIA_spde, theta = c(0,0))
  GIA_sprior <- sqrt(1/diag(Q00))
  proj<- inla.mesh.projector(mesh_GPS[[i]], projection = "longlat", dims = c(361,181))
  image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_sprior)), col = topo.colors(40),
             xlab = "Longitude", ylab = "Latitude", main = mains[i])
  points(GPSX, GPSY, cex = 0.5)
}
dev.off()

# ## Plot the results 3D
#vals <- GIA_sprior
#open3d()
#par3d(windowRect = c(100, 100, 900, 900))
#layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
#par3d(zoom = 0.8)
#t_lim <- range(vals)
#t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
#t_Cpal <- topo.colors(t_Clens, alpha=0)
#t_Cols<- t_Cpal[round((vals - t_lim[1])*100) + 1]
#plot3d(GPS_loc, pch = "20", size = 5, xlab = "", ylab = "", zlab = "", axe=FALSE)
#plot(mesh_GPS[[3]], rgl = TRUE, col= t_Cols, edge.color = rgb(0, 0.5, 0.6, alpha =0), add = TRUE)

## plot the color bar
#next3d(reuse = FALSE)
#bgplot3d({z=matrix(1:t_Clens, nrow = t_Clens)
#y=1
#x=seq(t_lim[1],t_lim[2],len = t_Clens)
#par(cex = 1.5, fin = c(8, 1), mai = c(0,0, 0.5, 0), oma = c(1, 0, 0, 0))
#image(x,y,z,col = topo.colors(t_Clens),axes=FALSE,xlab="",ylab="")
#title("INLA Posterior ICE6G GIA, mu = 0")
#axis(1)})

