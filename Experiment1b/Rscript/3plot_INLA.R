###################################################################
##                  Plot the INLA results                        ##
## This script produce plots of the hyperparameters and          ##
## the predicted GIA and uncertainties.                          ##
## The prediction map is done in ggplot2                         ##
###################################################################

library(ggplot2)
library(grid)
library(gridExtra)

#### 1 Plot hyperparameters
pars_GIA <- inla.spde2.result(res_inla, "GIA", GIA_spde, do.transf=TRUE)
theta_mean <- pars_GIA$summary.theta$mean
theta_sd <- pars_GIA$summary.theta$sd

## Find the mode of rho and sigma^2
lrho_mode <- pars_GIA$summary.log.range.nominal$mode
lrho_mean <- pars_GIA$summary.log.range.nominal$mean
lrho_sd <- pars_GIA$summary.log.range.nominal$sd
rho_mode <- exp(lrho_mean - lrho_sd^2)

lsigma_mode <- pars_GIA$summary.log.variance.nominal$mode
lsigma_mean <- pars_GIA$summary.log.variance.nominal$mean
lsigma_sd <- pars_GIA$summary.log.variance.nominal$sd
sigma_mode <- exp(lsigma_mean - lsigma_sd^2)

png(paste0(outdir, outname, "_hyperpar.png"), width = 960, height = 480, pointsize = 20)
par(mfrow = c(1,2))
plot(pars_GIA$marginals.range.nominal[[1]], type = "l",
     main = bquote(bold(rho("mode") == .(round(rho_mode, 4))))) # The posterior from inla output
plot(pars_GIA$marginals.variance.nominal[[1]], type = "l", xlim = c(0, 10),
     main = bquote(bold({sigma^2}("mode") == .(round(sigma_mode, 4))))) # The posterior from inla output

dev.off()



#### 2 Plot the predictions
### Assemble the dataset for plot
INLA_pred <- res_inla$summary.linear.predictor

### Get the predicted mean and uncertainties
pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
GPS_idx <- pred_idx[1:nrow(GPS_obs)]
GIA_idx <- pred_idx[-c(1:nrow(GPS_obs))]

## GPS 
GPS_u <- INLA_pred$sd[GPS_idx]
GPS_pred <- data.frame(lon = GPSx, lat = GPS_obs$lat, u = GPS_u)

## GIA
GIA_m <- INLA_pred$mean[GIA_idx] + GIA_mu
GIA_u <- INLA_pred$sd[GIA_idx]
proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(360,180), xlim = c(0, 359), ylim = c(-89.5, 89.5))
GIA_grid <- expand.grid(proj$x, proj$y)
GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                       mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                       u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))))
write.table(GIA_pred, file = paste0(outdir, outname, "_predict.txt"), row.names = FALSE, eol = "\r\n")

## Now plot
world_map <- map_data("world2")
## 1 The Prior
map_prior <- ggplot(data=GIA_prior) + geom_raster(aes(x = x_center, y = y_center, fill = trend)) + 
  coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(limits=c(0,359),  expand = c(0, 0)) + scale_y_continuous(limits=c(-89.5,89.5),  expand = c(0, 0)) + 
  scale_fill_gradientn(colors = terrain.colors(12), name = "mm/yr", limit = c(-15, 15),
                       guide = guide_colorbar(barwidth = 2, barheight = 10, label.position = "right", title.position = "bottom")) 
beauty <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = 'white'), 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        panel.border = element_blank())

map_prior <- map_prior + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                       colour="grey", fill = NA, alpha = 0.5) + beauty + ggtitle("Prior GIA mean field")

## 2 The prediction map
map_GIA <- ggplot(data=GIA_pred) + geom_raster(aes(x = lon, y = lat, fill = mean)) + 
  coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
  scale_fill_gradientn(colors = terrain.colors(12), name = "mm/yr", limit = c(-15, 15),
                       guide = guide_colorbar(barwidth = 2, barheight = 10, label.position = "right", title.position = "bottom")) 


map_GIA2 <- map_GIA + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                   colour="grey", fill = NA, alpha = 0.5)

map_GIA3 <- map_GIA2 + geom_point(data=GPS_pred, aes(x=lon, y=lat), pch=19, size = GPS_pred$u*3,
                                    col = "red", fill = "red", alpha=0.7) 

GIA_std <- subset(GIA_pred, lon %in% seq(0, 359, 15))
GIA_std <- subset(GIA_std, lat %in% seq(-82.5, 82.5, 10))

map_GIA4 <- map_GIA3 + geom_point(data=GIA_std, aes(x=lon, y=lat), pch=19, size = GIA_std$u*3,
                                  col = "grey", fill = "grey", alpha=0.5)

map_GIAf1 <- map_GIA2 + beauty + ggtitle("Predicted GIA") +
  scale_x_continuous(limits=c(0,360),  expand = c(0, 0)) + 
  scale_y_continuous(limits=c(-90,90),  expand = c(0, 0)) 

map_GIAf2 <- map_GIA4 + beauty + ggtitle("Predicted GIA and uncertaities") +
  scale_x_continuous(limits=c(0,360),  expand = c(0, 0)) + 
  scale_y_continuous(limits=c(-90,90),  expand = c(0, 0)) 

## 3 The uncertainty map 
map_unc <- ggplot(data=GIA_pred) + geom_raster(aes(x = lon, y = lat, fill = u)) + 
  coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
  scale_fill_gradientn(colors = terrain.colors(12), name = "mm/yr",
                       guide = guide_colorbar(barwidth = 2, barheight = 10, label.position = "right", title.position = "bottom")) 

map_unc2 <- map_unc + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                   colour="grey", fill = NA, alpha = 0.5)

map_unc3 <- map_unc2 + geom_point(data=GPS_pred, aes(x=lon, y=lat), pch=19, size = 0.5,
                                  col = "black", fill = "grey", alpha=0.3) 

map_uncf <- map_unc3 + beauty + ggtitle("Predicted uncertainties") +
  scale_x_continuous(limits=c(0,360),  expand = c(0, 0)) + 
  scale_y_continuous(limits=c(-90,90),  expand = c(0, 0)) 


## Zoom in for inspecting the uncertainties and correlation length.
## Choose the area near Alaska.
zoom_data <- subset(GIA_pred, lon > 180 & lon < 230 & lat > 50 & lat < 75)
map_zoom <- map_GIA3 +
  scale_x_continuous(limits=c(180,230),  expand = c(0.03, 0.03)) + scale_y_continuous(limits=c(50,75),  expand = c(0, 0)) 
## Calculate the approximate correlation length
cor_deg <- round(rho_mode /(2*pi) * 360)
cor_dist <- round(rho_mode * 6371)
map_zoom <- map_zoom +ggtitle("Zoom-in near Alaska", 
                              subtitle = paste("The estimated correlation length is about", eval(cor_dist), "km", "(", eval(cor_deg),"degree)"))

zoom_std <- subset(zoom_data, lon %in% seq(182,228, 2))
zoom_std <- subset(zoom_std, lat %in% seq(52.5, 72.5, 2))
map_zoom <- map_zoom + geom_point(data=zoom_std, aes(x=lon, y=lat), pch=19, size = zoom_std$u*3,
                                  col = "black", fill = "black", alpha=0.3) + beauty

## Plot on the grid  
png(paste0(outdir, outname, "_GIAmap1.png"), width = 900, height = 1200, pointsize = 20)
grid.arrange(map_prior, map_GIAf1, map_uncf)
dev.off()

png(paste0(outdir, outname, "_GIAmap2.png"), width = 900, height = 800, pointsize = 25)
grid.arrange(map_GIAf2, map_zoom)
dev.off()
