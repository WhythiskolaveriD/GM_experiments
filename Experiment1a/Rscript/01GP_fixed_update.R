###################################################################
##                    Gaussian Update                            ##
## This script performs a Gaussian Update of the GIA processes   ##
## Given                                                         ##
## 1. fixed parameter values                                     ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model prior                                    ##
###################################################################
library(rgl)
library(GEOmap)
library(INLA)
library(MVST)
library(fields)
library(gstat)
library(geoR)
load("C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/Mesh_GIA.RData")
GPS_obs <- read.table("Z:/ExperimentsBHM/Experiment1a/inputs/GPS_synthetic.txt", header = T)
GPS_obsU <- GPS_obs[!duplicated(GPS_obs[,2:3]), ]
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$y_center, lon = GPS_obs$x_center))

## Screening of data
## plot and calculate the sample variogram from longlat coords
coordinates(GPS_obsU) <- c("x_center", "y_center")
proj4string(GPS_obsU) <- CRS("+proj=longlat")
## calculate the great circle distance from long lat so set projection to be FALSE
GPS_obsU$trendc <- GPS_obsU$trend - fitted(lm(GPS_obsU$trend~x_center + y_center, data = GPS_obsU))
v1 <- variogram(trendc ~ 1, data = GPS_obsU)
plot(v1)
## From the plot, sill = sigma = 5, range = 500. choose kappa = 1 nugget = 0 (assume)
f1 <- fit.variogram(v1, vgm(5, model = "Mat", kappa = 1, range = 500))
summary(f1)  ## fitted range = 134.9, convert to unit ball assumming earth radius 6371km = 0.0212
plot(v1, f1)  ## Not a good fit from the plot


coordinates(GPS_obsU) <- c("x_center", "y_center")
proj4string(GPS_obsU) <- CRS("+proj=longlat")
GPSgeo <- as.geodata(GPS_obsU)
vv1 <- likfit(GPSgeo, ini.cov.pars = c(4, 10), kappa = 1, fix.nugget = TRUE, nugget = 0, 
              cov.model = "matern", lik.method = "REML")


## Better to use a weighted average of the prior info and estimates from observation -- Bayesian!

## Find the mapping between observations and processes basis
CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)

#### Fixed parameters for the precision matrix
## 1 Use Priors only
GIA_GMRF<- GMRF(mu = GIA_mu, Q = GIA_Q0) 
GPS_obsDF <- Obs(df=data.frame(x = GPS_loc[,1], y = GPS_loc[,2], w = GPS_loc[,3], z = GPS_obs$trend, 
                          std = GPS_obs$std, t = rep(0, nrow(GPS_loc))))
L <- link(GIA_GMRF, GPS_obsDF, Cmat = CMat)
e <- new("link_list",list(L))
v <- new("block_list",list(G1 = GIA_GMRF, O = GPS_obsDF))
G <- new("Graph",e = e,v = v)
G_reduced <- compress(G)
Results <- Infer(G_reduced)

GIA_post <- Results$Post_GMRF@rep
GIA_mpost<- GIA_post$x_mean
## Get the posterior mean and variance
GIA_spost <- sqrt(GIA_post$x_margvar)
GIA_diff <- GIA_mpost - GIA_mu

## 2 Use obs estimates only
range1 <- 0.0212
sigma1 <- sqrt(4.52)
kappa1 <- sqrt(8)/range0
tau1 <- 1/(4*pi*kappa1^2*sigma1^2)
GIA_Q1 <- inla.spde.precision(GIA_spde, theta=c(log(sqrt(tau1)), log(kappa1)))

GIA_GMRF1<- GMRF(Q = GIA_Q1) 
L1 <- link(GIA_GMRF1, GPS_obsDF, Cmat = CMat)
e1 <- new("link_list",list(L1))
v1 <- new("block_list",list(G1 = GIA_GMRF1, O = GPS_obsDF))
G1 <- new("Graph",e = e1,v = v1)
G_reduced1 <- compress(G1)
Results1 <- Infer(G_reduced1)

GIA_post1 <- Results1$Post_GMRF@rep
GIA_mpost1<- GIA_post1$x_mean
## Get the posterior mean and variance
GIA_spost1 <- sqrt(GIA_post1$x_margvar)
GIA_diff1 <- GIA_mpost1 - GIA_mu


## Plot the results
vals <- GIA_spost1
open3d()
par3d(windowRect = c(100, 100, 900, 900))
layout3d(matrix(1:2, 2,1), heights = c(4,2), sharedMouse = TRUE)
par3d(zoom = 0.8)
t_lim <- range(vals)
t_Clens <- round((t_lim[2] - t_lim[1])*100) + 1
t_Cpal <- topo.colors(t_Clens, alpha=0)
t_Cols<- t_Cpal[round((vals - t_lim[1])*100) + 1]
plot3d(GPS_loc, cex = 2, xlab = "", ylab = "", zlab = "", axe=FALSE)
plot(Mesh_GIA, rgl = TRUE, col= t_Cols, edge.color = rgb(0, 0.5, 0.6, alpha =0), add = TRUE)

## plot the color bar
next3d(reuse = FALSE)
bgplot3d({z=matrix(1:t_Clens, nrow = t_Clens)
y=1
x=seq(t_lim[1],t_lim[2],len = t_Clens)
par(cex = 1.5, fin = c(8, 1), mai = c(0,0, 0.5, 0), oma = c(1, 0, 0, 0))
image(x,y,z,col = topo.colors(t_Clens),axes=FALSE,xlab="",ylab="")
title("Posterior ICE6G GIA variance -- mu = 0")
axis(1)})

writeWebGL(dir = "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a", 
           filename= "C:/Users/zs16444/Local Documents/GlobalMass/Experiment1a/GIA_globe.html",  # save the plot as an html
           width = 1000, reuse = TRUE)

