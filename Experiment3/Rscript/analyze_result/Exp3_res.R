###### This script run the inla procedure for Experiment 3
####################################################################################################

#res_inla <- readRDS(file = "Z:/WP1-BHM/Experiment3/Exp3_inla.rds")

mass_sig <- inla.tmarginal(function(x) exp(x + lsig_massm), res_inla$marginals.hyperpar$`Theta1 for mass`)
mass_rho <- inla.tmarginal(function(x) exp(x + lrho_massm)*6371, res_inla$marginals.hyperpar$`Theta2 for mass`)

gia_sig <- inla.tmarginal(function(x) exp(x + lsig_giam), res_inla$marginals.hyperpar$`Theta1 for GIA`)
gia_rho <- inla.tmarginal(function(x) exp(x + lrho_giam)*6371, res_inla$marginals.hyperpar$`Theta2 for GIA`)

par(mfrow= c(2,2))
plot(mass_sig, type = "l", main = "marginal variance for mass (mm)")
plot(mass_rho, type = "l", main = "correlation length for mass (km)")
plot(gia_sig, type = "l", main = "marginal variance for gia (mm)")
plot(gia_rho, type = "l", main = "correlation length for gia (km)")

#### Assemble the results
## Extract the predictions
INLA_mass <- res_inla$summary.random$mass
INLA_GIA <- res_inla$summary.random$GIA

## mass
mass_mean <- INLA_mass$mean
mass_u <- INLA_mass$sd
proj <- inla.mesh.projector(mesh0, projection = "longlat", xlim = c(172, 328), ylim = c(-3, 89))
mass_grid <- expand.grid(proj$x, proj$y)
mass_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                        mean = as.vector(inla.mesh.project(proj, as.vector(mass_mean))),
                        u = as.vector(inla.mesh.project(proj, as.vector(mass_u))))

## gia
GIA_mean <- INLA_GIA$mean + meshgia
GIA_mean0 <- INLA_GIA$mean
GIA_u <- INLA_GIA$sd
GIA_pred <- data.frame(lon = mass_grid[,1], lat = mass_grid[,2],
                       mean = as.vector(inla.mesh.project(proj, as.vector(GIA_mean))),
                       mean0 = as.vector(inla.mesh.project(proj, as.vector(GIA_mean0))),
                       u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))))

ress_3 <- list(spde = list(mass = mass_spde, gia = gia_spde),  
               mesh = mesh0,
               pred = list(mass = mass_pred,gia = GIA_pred))


save(ress_3, file = paste0(dd_save, filename, "_res.RData"))

### Save ncdf files

## save mass
df2ncdf(df = mass_pred[, c("lon", "lat", "mean", "u")], 
        fname = paste0(dd_save, filename, "_mass.nc"),  
        vars = list(c("mean", "trend", "BHM predicted mass trend", "mm/yr"), 
                    c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
        title = "BHM-predicted mass", 
        append = FALSE)

## save GIA
df2ncdf(df = GIA_pred[,c("lon", "lat", "mean","u")],
        fname = paste0(dd_save,  filename, "_GIA.nc"),  
        vars = list(c("mean","trend", "BHM predicted GIA trend", "mm/yr"), 
                    c("u", "uncertainty", "uncertainty of the BHM prediction", "mm/yr")),
        title = "BHM-predicted GIA", 
        append = FALSE)

