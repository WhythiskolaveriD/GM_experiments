###################################################################
##          A Bayesian Gaussian Update (INLA)                    ##
## This script performs a full Bayesian Gaussian Update          ##
## of the GIA processes                                          ##
## Use a lot of memories -- run on nuuk at least                 ##
## Given                                                         ##
## 1. prior assumptions of the parameters                        ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model mean adjustment                          ##
###################################################################
#### INLA Approximation
## require
library(GEOmap)
library(INLA)
if(runserver){
INLA:::inla.dynload.workaround()
  }#Use this on server with old gcc library
library(fields)
library(maps)
library(maptools)
source("glbm/Experiment1a/Rscript/MVSTplus.R")

## Get the GPS locations
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$y_center, lon = GPS_obs$x_center))
GPSX <- ifelse(GPS_obs$x_center > 180, GPS_obs$x_center-360, GPS_obs$x_center)
GPSY <- GPS_obs$y_center
## Generate the mesh
mesh_GIA <- mesh_temp
MlocLL <- Lxyz2ll(list(x=mesh_GIA$loc[,1], y = mesh_GIA$loc[,2], z = mesh_GIA$loc[,3]))
MlocLL$lon <- ifelse(MlocLL$lon < 0, MlocLL$lon + 360, MlocLL$lon)
MlocLL$lon <- ifelse(MlocLL$lon > 359.5, MlocLL$lon - 360, MlocLL$lon)
M_sp <- SpatialPoints(data.frame(lon = MlocLL$lon, lat = MlocLL$lat), proj4string = CRS("+proj=longlat")) #This convert GIA_ice6g a SpatialPointDataFrame
Midx <- over(M_sp, Plist)
GIA_mu <- GIA_ice6g_sp$trend[Midx]


## Setup for spde parameters
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]

lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

GIA_spde <- inla.spde2.matern(mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))


## Find the mapping between observations and processes basis
Ay <- inla.spde.make.A(mesh = mesh_GIA, loc = GPS_loc)

## Create a stack of locations to plot the posterior (resolution = 5 degree apart)
s_lat <- seq(-89.5, 89.5, 10)
s_lon <- seq(-179.5, 179.5, 10)
sll <- expand.grid(x = s_lon, y = s_lat)
sll_loc <- do.call(cbind, Lll2xyz(lat = sll$y, lon = sll$x))
A_all <- inla.spde.make.A(mesh = mesh_GIA, loc = rbind(GPS_loc, mesh_GIA$loc, sll_loc))

## create new synthetic data, plot and compare
y0 <- as.vector(Ay %*% GIA_mu)
ydata <- GPS_obs$trend - y0

#### 1: GIA_ice6g prior mean info
st.est <- inla.stack(data = list(y=ydata), A = list(Ay),
                     effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
## Predict at the GPS location, mesh grid and predict error location
Ip <- Matrix(0, GIA_spde$n.spde, GIA_spde$n.spde)
diag(Ip) <- 1
st.pred <- inla.stack(data = list(y=NA), A = list(A_all),
                      effects = list(GIA=1:spde$n.spde), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

## Fix the sigma_e^2 to be 1 and scale them according to the data in inla(...,scale = scale)
## Default uses log(1/sigma_e^2) to be loggamma distribution with intial value = 0.
hyper <- list(prec = list(fixed = TRUE, initial = log(1/mean(err^2))))
formula = y ~ -1 + f(GIA, model = spde)
prec_scale <- 1/err^2
res_inla <- inla(formula, data = inla.stack.data(st.est, spde = spde), family = "gaussian",
                  control.family = list(hyper = hyper),
                   control.predictor=list(A=inla.stack.A(st.est), compute =TRUE))
summary(res_inla)

## Plot hyperparameters
pars_GIA <- inla.spde2.result(res_inla, "GIA", GIA_spde, do.transf=TRUE)
theta_mean <- pars_GIA$summary.theta$mean
theta_sd <- pars_GIA$summary.theta$sd

pdf(file = paste0(wkdir, exname, "hyperpar.pdf"), width = 6, height = 8)
par(mfrow = c(2,2))

## plot the range parameter rho
lrho_mode <- pars_GIA$summary.log.range.nominal$mode
lrho_mean <- pars_GIA$summary.log.range.nominal$mean
lrho_sd <- pars_GIA$summary.log.range.nominal$sd
rho_mode <- exp(lrho_mean - lrho_sd^2)

## plot log(rho)
plot(pars_GIA$marginals.log.range.nominal[[1]], type = "l",
     main = bquote(bold(log(rho)("mode") == .(round(lrho_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.log.range.nominal[[1]][,1]
yy0 <- dnorm(xx, mean = lrho0, sd = theta2_s) # The prior
lines(x = xx, y = yy0, col = 2)

## plot rho
plot(pars_GIA$marginals.range.nominal[[1]], type = "l",
     main = bquote(bold(rho("mode") == .(round(rho_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.range.nominal[[1]][,1]
yy0 <- dlnorm(xx, meanlog = lrho0, sdlog = sqrt(theta2_s)) # The prior
lines(x = xx, y = yy0, col = 2)


## Plot the field variance parameter sigma
lsigma_mode <- pars_GIA$summary.log.variance.nominal$mode
lsigma_mean <- pars_GIA$summary.log.variance.nominal$mean
lsigma_sd <- pars_GIA$summary.log.variance.nominal$sd
sigma_mode <- exp(lsigma_mean - lsigma_sd^2)

## plot log(sigma)
plot(pars_GIA$marginals.log.variance.nominal[[1]], type = "l",
     main = bquote(bold(log(sigma^2)("mode") == .(round(lsigma_mode, 4))))) # The posterior from inla output
xx <- pars_GIA$marginals.log.variance.nominal[[1]][,1]
yy0 <- dnorm(xx, mean = lsigma0, sd = theta1_s) # The prior
lines(x = xx, y = yy0, col = 2)

## plot sigma
plot(pars_GIA$marginals.variance.nominal[[1]], type = "l", xlim = c(0, 20),
     main = bquote(bold({sigma^2}("mode") == .(round(sigma_mode, 4))))) # The posterior from inla output
xx <- seq(from = 0.2, to=20, length.out = 1000)
yy0 <- dlnorm(xx, meanlog = lsigma0, sdlog = sqrt(theta1_s)) # The prior
lines(x = xx, y = yy0, col = 2)

dev.off()




## Plot the predicted GIA field mean and variance
proj <- inla.mesh.projector(mesh_GIA, projection = "longlat", dims = c(500,500))

pidx <- inla.stack.index(stGIA, tag = "pred")
n.data <- length(ydata)
n.GIA <- GIA_spde$n.spde
n.errs <- nrow(sll_loc)
GIA_mpost <- res_inla$summary.linear.predictor$mean[pidx$data[(n.data+1): (n.data+n.GIA)]] + GIA_mu
GIA_spost <- res_inla$summary.linear.predictor$sd[pidx$data[(n.data+1): (n.data+n.GIA)]]

GPS_mpost <- res_inla$summary.linear.predictor$mean[pidx$data[1:n.data]]
GPS_spost <- res_inla$summary.linear.predictor$sd[pidx$data[1:n.data]]

err_spost <- res_inla$summary.linear.predictor$sd[pidx$data[-c(1:(n.data+n.GIA))]]

GPScol <- ifelse(ydata > 0, 2, 1)

## plot the posterior error only
pdf(file = paste0(wkdir, exname, "GIA_error_field.pdf"), width = 8, height = 10)
par(mfrow = c(2,1))
Q2 <- inla.spde2.precision(GIA_spde, theta = theta_mean)
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(sqrt(1/diag(Q2)))), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field")
points(GPSX, GPSY, cex = GPS_spost*2, pch = 1)

## The standard error
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_spost)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Updated posterior error field")
points(GPSX, GPSY, cex = GPS_spost*2, pch = 1)
dev.off()


## plot the posterior error overlaid on posterior mean field
pdf(file = paste0(wkdir, exname, "GIAfield.pdf"), width = 15, height = 12)
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_mpost)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field")
points(sll$x, sll$y, cex = (err_spost/mean(err_spost))^2*2)
points(GPSX, GPSY, cex = (GPS_spost/mean(err_spost))^2*2, col = 2)
dev.off()


### plot zoom in
#indx <- which(proj$x > -120 & proj$x < -20)
#indy <- which(proj$y >0 & proj$y <70)
#par(mfrow = c(1,1))
#image.plot(proj$x[indx], proj$y[indy], inla.mesh.project(proj, as.vector(GIA_spost))[indx, indy], col = topo.colors(40),
#           xlab = "Longitude", ylab = "Latitude", main = "Updated posterior error field")
#points(sll$x, sll$y, cex = (err_spost/max(err_spost))^3*10)
#points(GPSX, GPSY, cex = GPS_spost^2*1.5, col = 2)


#image.plot(s_lon, s_lat, matrix(err_spost, nrow = length(s_lon), ncol = length(s_lat)), col = topo.colors(40),
#           xlab = "Longitude", ylab = "Latitude", main = "Updated posterior error field")
#points(sll$x, sll$y, cex = err_spost*2)
#points(GPSX, GPSY, cex = (GPS.sd/mean(GPS.sd))^2, col = 2)

save(res_inla, GIA_spde, ydata, mu_r, v_r, mu_s, v_s, mesh_GIA, file = paste0(wkdir, exname, "inla.RData"))



## Plot the results 3D
#vals <- GIA_spost
#open3d()
#par3d(windowRect = c(100, 100, 900, 900))
#layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
#par3d(zoom = 0.8)
#t_lim <- range(c(vals, err_spost))
#t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
#t_Cpal <- topo.colors(t_Clens, alpha=0)
#t_Cols<- t_Cpal[round((vals - t_lim[1])*100) + 1]
#GPScol <- t_Cpal[round((err_spost - t_lim[1])*100) + 1]
#plot3d(sll_loc, pch = "20", size = 5, col = GPScol, xlab = "", ylab = "", zlab = "", axe=FALSE)
#plot(mesh_GIA, rgl = TRUE, col= t_Cols, edge.color = rgb(0, 0.5, 0.6, alpha =0), add = TRUE)

## plot the color bar
#next3d(reuse = FALSE)
#bgplot3d({z=matrix(1:t_Clens, nrow = t_Clens)
#y=1
#x=seq(t_lim[1],t_lim[2],len = t_Clens)
#par(cex = 1.5, fin = c(8, 1), mai = c(0,0, 0.5, 0), oma = c(1, 0, 0, 0))
#image(x,y,z,col = topo.colors(t_Clens),axes=FALSE,xlab="",ylab="")
#title("INLA Posterior ICE6G GIA, mu = 0")
#axis(1)})




