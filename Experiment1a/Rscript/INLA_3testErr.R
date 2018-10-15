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

## Get the GPS locations
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$y_center, lon = GPS_obs$x_center))
GPSX <- ifelse(GPS_obs$x_center > 180, GPS_obs$x_center-360, GPS_obs$x_center)
GPSY <- GPS_obs$y_center

## load mesh and build spde model
load("experimentBHM/GIA_sp_info.RData")
load("experimentBHM/mesh_reg.RData")
Mesh_GIA <- mesh_regs[[1]]
MlocLL <- Lxyz2ll(list(x=Mesh_GIA$loc[,1], y = Mesh_GIA$loc[,2], z = Mesh_GIA$loc[,3]))
MlocLL$lon <- ifelse(MlocLL$lon < 0, MlocLL$lon + 360, MlocLL$lon)
MlocLL$lon <- ifelse(MlocLL$lon > 359.5, MlocLL$lon - 360, MlocLL$lon)
M_sp <- SpatialPoints(data.frame(lon = MlocLL$lon, lat = MlocLL$lat), proj4string = CRS("+proj=longlat")) #This convert GIA_ice6g a SpatialPointDataFrame
Midx <- over(M_sp, Plist)
GIA_mu <- GIA_ice6g_sp$trend[Midx]

## Setup for spde parameters
source("glbm/Experiment1a/Rscript/MVSTplus.R")
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))



## Find the mapping between observations and processes basis
Ay <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)

## Create a stack of locations to plot the posterior (resolution = 5 degree apart)
s_lat <- seq(-89.5, 89.5, 10)
s_lon <- seq(-179.5, 179.5, 10)
sll <- expand.grid(x=s_lon, y = s_lat)
sll_loc <- do.call(cbind, Lll2xyz(lat = sll$y, lon = sll$x))
As <- inla.spde.make.A(mesh = Mesh_GIA, loc = sll_loc)

## create new synthetic data, plot and compare
y0 <- as.vector(Ay %*% GIA_mu)
ydata <- GPS_obs$trend - y0
err1 <- err2 <- GPS_obs$std
err1[459] <- err1[459]*0.01
err2[459] <- err2[459]*100


## Build the estimation stack
st.est <- inla.stack(data = list(y=ydata), A = list(Ay), 
                      effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")

## Predict at the GPS location, mesh grid and predict error location
Ip <- Matrix(0, GIA_spde$n.spde, GIA_spde$n.spde)
diag(Ip) <- 1
st.pred <- inla.stack(data = list(y=NA), A = list(rbind(Ay, Ip, As)), 
                      effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
stGIA <- inla.stack(st.est, st.pred)

## Fix the sigma_e^2 to be 1 and scale them according to the data in inla(...,scale = scale)
## Default uses log(1/sigma_e^2) to be loggamma distribution with intial value = 0.
hyper <- list(prec = list(fixed = TRUE, initial = 0))
formula = y ~ -1 + f(GIA, model = GIA_spde)
prec_scale1 <- c(1/err1, rep(1, length(ydata)), rep(1, GIA_spde$n.spde), rep(1, nrow(sll_loc)))
prec_scale2 <- c(1/err2, rep(1, length(ydata)), rep(1, GIA_spde$n.spde), rep(1, nrow(sll_loc)))
res_inla1 <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian", 
                 scale = prec_scale1, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))

res_inla2 <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian", 
                  scale = prec_scale2, control.family = list(hyper = hyper),
                  control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))
summary(res_inla1)
summary(res_inla2)


## Plot hyperparameters
pars_GIA1 <- inla.spde2.result(res_inla1, "GIA", GIA_spde, do.transf=TRUE)
pars_GIA2 <- inla.spde2.result(res_inla2, "GIA", GIA_spde, do.transf=TRUE)

pdf(file = paste0(wkdir, exname, "hyperpar.pdf"), width = 6, height = 8)
par(mfrow = c(2,2))
## plot log(rho)
plot(pars_GIA1$marginals.log.range.nominal[[1]], type = "l", 
     main =expression(bold(log(rho)))) # The posterior from inla output
lines(pars_GIA2$marginals.log.range.nominal[[1]], type = "l", col = 2) # The posterior from inla output

## plot rho
plot(pars_GIA1$marginals.range.nominal[[1]], type = "l", 
     main = expression(bold(rho))) # The posterior from inla output
lines(pars_GIA2$marginals.range.nominal[[1]], type = "l", col = 2) 

## plot log(sigma)
plot(pars_GIA1$marginals.log.variance.nominal[[1]], type = "l",
     main = expression(bold(log(sigma^2)))) # The posterior from inla output
lines(pars_GIA2$marginals.log.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output

## plot sigma
plot(pars_GIA1$marginals.variance.nominal[[1]], type = "l", xlim = c(0, 20), 
     main = expression(bold(sigma^2))) # The posterior from inla output
lines(pars_GIA2$marginals.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output

dev.off()


## Plot the predicted GIA field mean and variance
proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(500,500))

pidx <- inla.stack.index(stGIA, tag = "pred")
n.data <- length(ydata)
n.GIA <- GIA_spde$n.spde
n.errs <- nrow(sll_loc)
GIA_mpost1 <- res_inla1$summary.linear.predictor$mean[pidx$data[(n.data+1): (n.data+n.GIA)]] + GIA_mu
GPS_spost1 <- res_inla1$summary.linear.predictor$sd[pidx$data[1:n.data]]
err_spost1 <- res_inla1$summary.linear.predictor$sd[pidx$data[-c(1:(n.data+n.GIA))]]

GIA_mpost2 <- res_inla2$summary.linear.predictor$mean[pidx$data[(n.data+1): (n.data+n.GIA)]] + GIA_mu
GPS_spost2 <- res_inla2$summary.linear.predictor$sd[pidx$data[1:n.data]]
err_spost2 <- res_inla2$summary.linear.predictor$sd[pidx$data[-c(1:(n.data+n.GIA))]]

## plot the posterior error overlaid on posterior mean field
pdf(file = paste0(wkdir, exname, "GIAfield.pdf"), width = 12, height = 18)
par(mfrow= c(2,1))
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_mpost1)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field")
points(sll$x, sll$y, cex = (err_spost1/mean(err_spost1))^2*2)
points(GPSX, GPSY, cex = (GPS_spost1/mean(err_spost1))^2*2, col = 2)
rect(xleft =-15 , ybottom =-85 , xright = 10, ytop = -60, lwd = 1.5, border = 2)

image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_mpost2)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field")
points(sll$x, sll$y, cex = (err_spost2/mean(err_spost2))^2*2)
points(GPSX, GPSY, cex = (GPS_spost2/mean(err_spost2))^2*2, col = 2)
rect(xleft =-15 , ybottom =-85 , xright = 10, ytop = -60, lwd = 1.5, border = 2)

dev.off()



save(res_inla1, res_inla2, GIA_spde, ydata, err1, err2, mu_r, v_r, mu_s, v_s, Mesh_GIA, file = paste0(wkdir, exname, "inla.RData"))


