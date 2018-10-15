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
## Create store objs
res_inlas <- list()
projs <- list()
pars_GIAs <- list()
GIA_mposts <- list()
err_sposts <- list()
GPS_sposts <- list()

## Define mesh size
meshSize <- c("s", "m", "l")
## load GPS data
load("experimentBHM/GIA_sp_info.RData")
GPS_obs <- read.table("experimentBHM/GPS_synthetic.txt", header = T)
load("experimentBHM/mesh_reg.RData")

## use the regular mesh
for (i in 1:3){
  exname <- paste0("/experimentBHM/reg_",meshSize[i],"Mesh")
  mesh_temp <- mesh_regs[[i]]
  source("glbm/Experiment1a/Rscript/INLA_Bayes_update.R")
  res_inlas[[i]] <- res_inla
  pars_GIAs[[i]] <- pars_GIA
  projs[[i]] <- proj
  GIA_mposts[[i]] <- GIA_mpost
  err_sposts[[i]] <- err_spost
  GPS_sposts[[i]] <- GPS_spost
}


## Use the irregular mesh
load("experimentBHM/mesh_GPS.RData")
exname <- paste0("/experimentBHM/irreg_",meshSize[i],"Mesh")
mesh_temp <- mesh_GPS[[3]]
source("glbm/Experiment1a/Rscript/INLA_Bayes_update.R")
res_inlas[[4]] <- res_inla
pars_GIAs[[4]] <- pars_GIA
projs[[4]] <- proj
GIA_mposts[[4]] <- GIA_mpost
err_sposts[[4]] <- err_spost
GPS_sposts[[4]] <- GPS_spost

## Plot hyperparameters -- small errors
exname <- paste0("/experimentBHM/","Meshcompare")
pdf(file = paste0(wkdir, exname, "_hyperpar.pdf"), width = 6, height = 8)
par(mfrow = c(2,2))
## plot log(rho)
maxd <- max(c(pars_GIAs[[1]]$marginals.log.range.nominal[[1]][,2],
              pars_GIAs[[2]]$marginals.log.range.nominal[[1]][,2],
              pars_GIAs[[3]]$marginals.log.range.nominal[[1]][,2],
              pars_GIAs[[4]]$marginals.log.range.nominal[[1]][,2]))
plot(pars_GIAs[[1]]$marginals.log.range.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main =expression(bold(log(rho)))) # The posterior from inla output
lines(pars_GIAs[[2]]$marginals.log.range.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3]]$marginals.log.range.nominal[[1]], type = "l", col = 4)
lines(pars_GIAs[[4]]$marginals.log.range.nominal[[1]], type = "l", col = 5)

## plot rho
maxd <- max(c(pars_GIAs[[1]]$marginals.range.nominal[[1]][,2],
              pars_GIAs[[2]]$marginals.range.nominal[[1]][,2],
              pars_GIAs[[3]]$marginals.range.nominal[[1]][,2],
              pars_GIAs[[4]]$marginals.range.nominal[[1]][,2]))
plot(pars_GIAs[[1]]$marginals.range.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main = expression(bold(rho))) # The posterior from inla output
lines(pars_GIAs[[2]]$marginals.range.nominal[[1]], type = "l", col = 2) 
lines(pars_GIAs[[3]]$marginals.range.nominal[[1]], type = "l", col = 4)
lines(pars_GIAs[[4]]$marginals.range.nominal[[1]], type = "l", col = 5)


## plot log(sigma)
maxd <- max(c(pars_GIAs[[1]]$marginals.log.variance.nominal[[1]][,2],
              pars_GIAs[[2]]$marginals.log.variance.nominal[[1]][,2],
              pars_GIAs[[3]]$marginals.log.variance.nominal[[1]][,2],
              pars_GIAs[[4]]$marginals.log.variance.nominal[[1]][,2]))

plot(pars_GIAs[[1]]$marginals.log.variance.nominal[[1]], type = "l", ylim = c(0,signif(maxd,2)),
     main = expression(bold(log(sigma)))) # The posterior from inla output
lines(pars_GIAs[[2]]$marginals.log.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3]]$marginals.log.variance.nominal[[1]], type = "l", col = 4) # The posterior from inla output
lines(pars_GIAs[[4]]$marginals.log.variance.nominal[[1]], type = "l", col = 5) # The posterior from inla output

## plot sigma
maxd <- max(c(pars_GIAs[[1]]$marginals.variance.nominal[[1]][,2],
              pars_GIAs[[2]]$marginals.variance.nominal[[1]][,2],
              pars_GIAs[[3]]$marginals.variance.nominal[[1]][,2],
              pars_GIAs[[4]]$marginals.variance.nominal[[1]][,2]))
maxv <- max(c(pars_GIAs[[1]]$marginals.variance.nominal[[1]][,1],
              pars_GIAs[[2]]$marginals.variance.nominal[[1]][,1],
              pars_GIAs[[3]]$marginals.variance.nominal[[1]][,1],
              pars_GIAs[[4]]$marginals.variance.nominal[[1]][,1]))

plot(pars_GIAs[[1]]$marginals.variance.nominal[[1]], type = "l",  xlim = c(0, maxv/20), ylim = c(0, signif(maxd, 2)),
     main = expression(bold(sigma))) # The posterior from inla output
lines(pars_GIAs[[2]]$marginals.variance.nominal[[1]], type = "l", col = 2) # The posterior from inla output
lines(pars_GIAs[[3]]$marginals.variance.nominal[[1]], type = "l", col = 4) # The posterior from inla output
lines(pars_GIAs[[4]]$marginals.variance.nominal[[1]], type = "l", col = 5) # The posterior from inla output
legend("topright", legend = c("reg_s", "reg_m", "reg_l", "GPS_l"), lwd = rep(1,4), lty = rep(1,4), col = c(1, 2, 4, 5))
dev.off()


## Plot the predicted GIA field mean and variance
pdf(file = paste0(wkdir, exname, "_GIAfield.pdf"), width = 10, height = 10)
par(mfrow = c(2,1))
image.plot(projs[[1]]$x, projs[[1]]$y, inla.mesh.project(projs[[1]], as.vector(GIA_mposts[[1]])), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field -- reg small")
points(sll$x, sll$y, cex = (err_sposts[[1]]/mean(err_sposts[[1]]))^2*2)
points(GPSX, GPSY, cex = (GPS_spost[[1]]/mean(err_sposts[[1]]))^2*2, col = 2)

image.plot(projs[[2]]$x, projs[[2]]$y, inla.mesh.project(projs[[2]], as.vector(GIA_mposts[[2]])), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field -- reg medium")
points(sll$x, sll$y, cex = (err_sposts[[2]]/mean(err_sposts[[2]]))^2*2)
points(GPSX, GPSY, cex = (GPS_spost[[2]]/mean(err_sposts[[2]]))^2*2, col = 2)

par(mfrow = c(2,1))
image.plot(projs[[2]]$x, projs[[3]]$y, inla.mesh.project(projs[[3]], as.vector(GIA_mposts[[3]])), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field -- reg large")
points(sll$x, sll$y, cex = (err_sposts[[3]]/mean(err_sposts[[3]]))^2*2)
points(GPSX, GPSY, cex = (GPS_spost[[3]]/mean(err_sposts[[3]]))^2*2, col = 2)

image.plot(projs[[4]]$x, projs[[4]]$y, inla.mesh.project(projs[[4]], as.vector(GIA_mposts[[4]])), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field -- irreg large")
points(sll$x, sll$y, cex = (err_sposts[[4]]/mean(err_sposts[[4]]))^2*2)
points(GPSX, GPSY, cex = (GPS_spost[[4]]/mean(err_sposts[[4]]))^2*2, col = 2)


dev.off()
