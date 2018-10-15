if(Sys.info()["nodename"] == "IT034844"){
  ## BS
  wd <- "C:/ZSwork/"
  dd <- "Z:/"
}else if(Sys.info()["nodename"] == "IT064613"){
  ## Maths
  wd <-"/home/zs16444/coding/"
  dd <-"/home/zs16444/globalmass/"
}else{
  ## server  
  wd <- "~/"
  dd <- "/./projects/GlobalMass/"
}

library(rgdal); library(sp);library(GEOmap)
library(INLA); library(ggplot2)

#### Simulate some data
## regular grids on [0,5]
loc_grid <- as.matrix(expand.grid(seq(0.05, 5, 0.1), seq(0.05, 5, 0.1)))

## Calculate the Euclidean distance matrix
locdist <- as.matrix(dist(loc_grid))

## Calculate the Matern covariance function for these grid points
Matern_cov <- fields::Matern(d =locdist, range = 0.5, smoothness = 1, phi = 2) 
## -- range = rho, the correlation length
## -- smoothness is nu as in the SPDE approach, nu = alpha - d/2, 
## -- alpha = 2 by default in INLA, d = 2 for 2D, so nu = 1 here
## -- phi is the marginal variance

## Simulate a realization of this process on the grid
set.seed(10)
data0 <- rnorm(nrow(locdist))
Matern_chol<- chol(Matern_cov) 
X <- Matern_chol %*% data0
Xdf <- data.frame(x = loc_grid[,1], y = loc_grid[,2], Xproc = X)

## Convert to a spatial grid data
Xsp <- Xdf
coordinates(Xsp) <- ~x+y
spplot(Xsp)

## Add some error 
err <- rnorm(nrow(loc_grid))*0.2
Xsp$obs <- Xsp$Xproc+err
spplot(Xsp)

## Use INLA
## Generate the mesh
mesh <- inla.mesh.2d(loc = as.matrix(expand.grid(seq(0, 5, 0.3), seq(0, 5, 0.3))), 
                     offset = c(0.6, 0.6), max.edge = c(0.5, 0.5))

## set prior values for the hyper parameters, delibrately set to be different from the truth.
sigma0 <- 1
range0 <- 1
lkappa0 <- log(8)/2 - log(range0)
ltau0 <- 0.5*log(1/(4*pi)) - log(sigma0) - lkappa0

## build the spde model
spde <- inla.spde2.matern(mesh, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3), 
                          theta.prior.mean = c(0,0), theta.prior.prec = c(0.01, 0.01))

## build the data stack
data_loc <- as.matrix(loc_grid)
data_obs <- Xsp$obs
Ay <- inla.spde.make.A(mesh = mesh, loc = data_loc)
st.est <- inla.stack(data = list(y=data_obs), A = list(Ay),
                     effects = list(process = 1:spde$n.spde), tag = "est")

hyper <- list(prec = list(fixed = TRUE, initial = log(1)))
formula = y ~ -1 + f(process, model = spde)
prec_scale <- (1/0.2)^2

## Generate 3 regions for aggregation
tri_sp <- SpatialPoints(coords = mesh$loc[,1:2])
region1 <- which(mesh$loc[,1] > 0.5 & mesh$loc[,1] < 3 & mesh$loc[,2] > 0.5 & mesh$loc[,2] < 2)
region2 <- which(mesh$loc[,1] > 3.2 & mesh$loc[,1] < 4.5 & mesh$loc[,2] > 0.5 & mesh$loc[,2] < 4)
region3 <- which(mesh$loc[,1] > 0.2 & mesh$loc[,1] < 3 & mesh$loc[,2] > 2.5 & mesh$loc[,2] < 4.5)
col <- rep(1, 578)
col[region1] <- 2
col[region2] <- 3
col[region3] <- 4
plot(tri_sp, col = col)
## predict a linear compbination 
A <- matrix(0, nrow = 3, ncol= 578)
A[1, region1] <- 1/length(region1)
A[2, region2] <- 1/length(region2)
A[3, region3] <- 1/length(region3)

lc.A = inla.make.lincombs(process=A)


res_inla <- inla(formula, data = inla.stack.data(st.est, spde = spde), family = "gaussian",
                 lincomb = lc.A,
                 scale = prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(st.est), compute =TRUE),
                 control.inla = list(lincomb.derived.only=FALSE), verbose=TRUE)


res_inla$summary.lincomb.derived$sd
sqrt(diag(A %*% Diagonal(x=res_inla$summary.random$process$sd^2) %*% t(A)))






## the predictions
proj <- inla.mesh.projector(mesh,dims = c(500,500))
mean_post <- res_inla$summary.random$process$mean
sd_post <- res_inla$summary.random$process$sd

xygrid <- expand.grid(proj$x, proj$y)
predmean <- as.vector(inla.mesh.project(proj, as.vector(mean_post)))
predu <- as.vector(inla.mesh.project(proj, as.vector(sd_post)))
predata <- data.frame(x = xygrid[,1], y = xygrid[,2], mean = predmean, u = predu)

ggplot(predata) + coord_fixed() + scale_x_continuous(limits=c(0,5),  expand = c(0.02, 0.02)) + 
  scale_y_continuous(limits=c(0,5),  expand = c(0.02, 0.02)) + geom_raster(aes(x = x, y = y, fill = mean))  + 
  scale_fill_gradientn(colours = terrain.colors(12), limit = c(-5,5))+ ggtitle("predicted mean")

ggplot(predata) + coord_fixed() + scale_x_continuous(limits=c(0,5),  expand = c(0.02, 0.02)) + 
  scale_y_continuous(limits=c(0,5),  expand = c(0.02, 0.02)) + geom_raster(aes(x = x, y = y, fill = u))  +
  scale_fill_gradientn(colours = terrain.colors(12), limit = c(0,0.3))+ ggtitle("predicted uncertainty")
