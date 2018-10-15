### This script show summary of the results and do some basic plots

## First load the saved results fi you started from a new session; otherwise ignore the following chunck
## Change the file names accordingly
dd_root <- "Z:/"
filename <- "Exp2_GSFCgrid_raw"
load(paste0(dd_root, "WP1-BHM/Experiment2/", filename, "_res.RData"))
load(paste0(dd_root, "WP1-BHM/Experiment2/data/Exp2_GSFCgrid_stks.rda"))
library(INLA)


## Summary the inference results
res_inla <- ress_2$res_inla
summary(res_inla)

## Plot the posteriors of the hyper parameters
## log(r) = log(r0) + theta2
## log(s) = log(s0) + theta1

## The posterior modes
## Mass
## ocean
ocean_sig <- inla.tmarginal(function(x) exp(x + lsig_oceanm), res_inla$marginals.hyperpar$`Theta1 for mass`)
ocean_rho <- inla.tmarginal(function(x) exp(x + lrho_oceanm)*6371, res_inla$marginals.hyperpar$`Theta2 for mass`)

## land
land_sig <- inla.tmarginal(function(x) exp(x + lsig_landm), res_inla$marginals.hyperpar$`Theta3 for mass`)
land_rho <- inla.tmarginal(function(x) exp(x + lrho_landm)*6371, res_inla$marginals.hyperpar$`Theta4 for mass`)

## Steric
steric_sig <- inla.tmarginal(function(x) exp(x + lsigma2), res_inla$marginals.hyperpar$`Theta1 for steric`)
steric_rho <- inla.tmarginal(function(x) exp(x + lrho2)*6371, res_inla$marginals.hyperpar$`Theta2 for steric`)

par(mfrow= c(2,2))
plot(ocean_sig, type = "l", main = "marginal variance for ocean mass (mm)")
plot(ocean_rho, type = "l", main = "correlation length for ocean mass (km)")
plot(land_sig, type = "l", main = "marginal variance for land mass (mm)")
plot(land_rho, type = "l", main = "correlation length for land mass (km)")

par(mfrow = c(1,2))
plot(steric_sig, type = "l", main = "marginal variance for steric (mm)")
plot(steric_rho, type = "l", main = "correlation length for steric (km)")



## Check that mass sum to zero as the constraints in the model run
## constraint 1 -- vertices sum to zero
sum(res_inla$summary.random$mass$mean)

## Constraint 2 -- grace loc sum to zero
sum(A_grace %*% res_inla$summary.random$mass$mean)

## Constraint 3 -- longlat grid sum to zero
sum(mass_pred$mean * mass_pred$areas)/sum(mass_pred$areas)


## Plot the predictions
steric_pred <- ress_2$pred$steric
ssh_pred <- ress_2$pred$ssh
mass_pred <- ress_2$pred$mass

## plot the mass mean 
lattice::levelplot(mean2 ~ lon + lat, data = mass_pred, aspect = "iso", at = seq(-20, 20, 2),
                   panel = function(x,y,z,...){
                     lattice::panel.fill(col = "grey")
                     lattice::panel.levelplot(x,y,z,...)
                     map2 <- map("world2", interior = FALSE, plot = FALSE)
                     lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                   },
                   main = "The predicted mass trend (mm/yr ewh)", xlab = "longitude", ylab = "latitude")

## plot the mass uncertainty
lattice::levelplot(u2 ~ lon + lat, data = mass_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
                   panel = function(x,y,z,...){
                     lattice::panel.fill(col = "grey")
                     lattice::panel.levelplot(x,y,z,...)
                     map2 <- map("world2", interior = FALSE, plot = FALSE)
                     lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                   },
                   main = "The predited mass uncertainties (mm/yr ewh)", xlab = "longitude", ylab = "latitude")

## Plot the steric mean
lattice::levelplot(mean2 ~ lon + lat, data = steric_pred, aspect = "iso", at = seq(-20, 20, 2),
                   panel = function(x,y,z,...){
                     lattice::panel.fill(col = "grey")
                     lattice::panel.levelplot(x,y,z,...)
                     map2 <- map("world2", interior = FALSE, plot = FALSE)
                     lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                   },
                   main = "The predicted steric trend (mm/yr)", xlab = "longitude", ylab = "latitude")

## plot the steric uncertainty
lattice::levelplot(u2 ~ lon + lat, data = steric_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
                   panel = function(x,y,z,...){
                     lattice::panel.fill(col = "grey")
                     lattice::panel.levelplot(x,y,z,...)
                     map2 <- map("world2", interior = FALSE, plot = FALSE)
                     lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                   },
                   main = "The predited steric uncertainties (mm/yr)", xlab = "longitude", ylab = "latitude")

## plot the ssh mean 
lattice::levelplot(mean2 ~ lon + lat, data = ssh_pred, aspect = "iso", at = seq(-20, 20, 1),
                   panel = function(x,y,z,...){
                     lattice::panel.fill(col = "grey")
                     lattice::panel.levelplot(x,y,z,...)
                     map2 <- map("world2", interior = FALSE, plot = FALSE)
                     lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                   },
                   main = "The predicted SSH trend (mm/yr)", xlab = "longitude", ylab = "latitude")

## plot the ssh uncertainty
lattice::levelplot(u2 ~ lon + lat, data = ssh_pred, aspect = "iso", at = seq(0, 4, 0.5),col.regions = topo.colors(10),
                   panel = function(x,y,z,...){
                     lattice::panel.fill(col = "grey")
                     lattice::panel.levelplot(x,y,z,...)
                     map2 <- map("world2", interior = FALSE, plot = FALSE)
                     lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                   },
                   main = "The predited SSH uncertainties (mm/yr)", xlab = "longitude", ylab = "latitude")



# Predict steric on Basin level 
load(paste0(dd, "WP1-BHM/maps/Ocean/basins.rda"))
pred_basin <- ress_2$pred$basins
pred_basin$names <- c("SAt", "ISPac", "EPac", "STNAt", "SPNAt", "NWPac", "global")

allbasinsdf$trend1 <- pred_basin$mean[1:6]
allbasinsdf$trend2 <- pred_basin$meanpost[1:6]
allbasinsdf$u1 <- pred_basin$u[1:6]
allbasinsdf$u2 <- pred_basin$u2[1:6]
allbasinsdf$u3 <- pred_basin$upost[1:6]

spplot(allbasinsdf, "trend1", at = c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5), col.regions = rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4')),
       main = "Basin trend (posterior mean)")

spplot(allbasinsdf, "trend2", at = c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5), col.regions = rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4')),
       main = "Basin trend (posterior sample mean)")

spplot(allbasinsdf, "u1", at = c(0, 0.04, 0.08, 0.12, 0.16, 0.2), 
       col.regions = c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000'),
       main = "uncertainty (posterior covariance approximation)")

spplot(allbasinsdf, "u2", at = c(0, 0.04, 0.08, 0.12, 0.16, 0.2), 
       col.regions = c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000'),
       main = "uncertainty (independent errors)")

spplot(allbasinsdf, "u3", at = c(0, 0.04, 0.08, 0.12, 0.16, 0.2), 
       col.regions = c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000'),
       main = "uncertainty (posterior sample approximation)")


