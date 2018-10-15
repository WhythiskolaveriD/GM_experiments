## summary and plot legacy

## print plots
my.dev.print <- function(path, filename, format = c("pdf", "pgn")){
  if(format == "png"){
    fname <- paste0(path, filename, ".png") 
    dev.print(png, fname, width = par("din")[1], height = par("din")[2], 
              units = "in", pointsize = 12, res = 72, bg = "transparent")
  }else if(format == "pdf"){
    fname <- paste0(path, filename, ".pdf") 
    dev.print(pdf, fname)
  }
  
}


## Wrapper for esitmation and prediction using INLA  
BayesDA_GIA<- function(Mesh_GIA, GPS_data){
  GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                                theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))
  
  ## Link the process to observations and predictions
  A_data <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
  A_pred <- inla.spde.make.A(mesh = Mesh_GIA, loc = rbind(GPS_loc, Mesh_GIA$loc))
  
  ## Create the estimation and prediction stack
  st.est <- inla.stack(data = list(y=GPS_data$trend0), A = list(A_data),
                       effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
  st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                        effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
  stGIA <- inla.stack(st.est, st.pred)
  
  ## Fix the GPS errors
  prec_scale <- c(1/GPS_data$std^2, rep(1, nrow(A_pred)))
  
  ## Run INLA
  res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian",
                   scale =prec_scale, control.family = list(hyper = hyper),
                   control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))
  
  INLA_pred <- res_inla$summary.linear.predictor
  
  ## Extract and project predictions
  pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
  GPS_idx <- pred_idx[1:nrow(GPS_data)]
  GIA_idx <- pred_idx[-c(1:nrow(GP_dataS))]
  
  ## GPS 
  GPS_u <- INLA_pred$sd[GPS_idx]
  GPS_pred <- data.frame(lon = GPSX, lat = GPS$lat, u = GPS_u)
  
  ## GIA
  GIA_diff <- INLA_pred$mean[GIA_idx] 
  GIA_m <- GIA_diff + GIA_prior
  GIA_u <- INLA_pred$sd[GIA_idx]
  proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(360,180), xlim = c(0, 360), ylim = c(-90, 90))
  GIA_grid <- expand.grid(proj$x, proj$y)
  GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                         diff = as.vector(inla.mesh.project(proj, as.vector(GIA_diff))),
                         mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                         u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))))
  
  return(list(res_inla = res_inla, spde = GIA_spde, st = stGIA, 
              mesh = Mesh_GIA, GPS_pred = GPS_pred, GIA_pred = GIA_pred))
}


BayesDA<- function(GIA, GPS, trho, tsigma){
  
  fibo_points <- fiboSphere2(N = 12960)
  mesh_points_xyz <- do.call(cbind, Lll2xyz(lat = fibo_points[,2], lon = fibo_points[,1]))
  Mesh_GIA <- inla.mesh.2d(loc = mesh_points_xyz, cutoff = 0.01, max.edge = 0.5)
  
  ## Transform the GIA data
  polycoords <- GIA[,c(6:13, 6,7)] 
  plist <- lapply(GIA$ID, function(x) Polygons(list(Polygon(cbind(lon = as.numeric(polycoords[x, c(1,3,5,7,9)]), 
                                                                  lat = as.numeric(polycoords[x, c(2,4,6,8,10)])))), ID = x))
  Plist <- SpatialPolygons(plist, proj4string = CRS("+proj=longlat"))
  ## Convert the Cartesian to longlat with long in [0, 360)
  MlocLL <- Lxyz2ll(list(x=Mesh_GIA$loc[,1], y = Mesh_GIA$loc[,2], z = Mesh_GIA$loc[,3]))
  MlocLL$lon <- ifelse(MlocLL$lon <0, MlocLL$lon + 359, MlocLL$lon)
  M_sp <- SpatialPoints(data.frame(lon = MlocLL$lon, lat = MlocLL$lat), proj4string = CRS("+proj=longlat")) 
  Midx <- over(M_sp, Plist)
  GIA_prior <- GIA$trend[Midx]
  
  ## Transform the GPS data
  GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS$lat, lon = GPS$lon))
  GPSX <- ifelse(GPS$lon < 0, GPS$lon + 359, GPS$lon)
  GPS_sp <- SpatialPoints(data.frame(lon = GPSX, lat = GPS$lat), proj4string = CRS("+proj=longlat"))
  Midx2 <- over(GPS_sp, Plist)
  GPS_mu <- GIA$trend[Midx2]
  GPS_data <- GPS$trend - GPS_mu
  
  ## The SPDE model
  lsigma0 <- tsigma[1]
  theta1_s <- tsigma[2]
  lrho0 <- trho[1]
  theta2_s <- trho[2]
  lkappa0 <- log(8)/2 - lrho0
  ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0
  GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                                theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))
  
  ## Link the process to observations and predictions
  A_data <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
  A_pred <- inla.spde.make.A(mesh = Mesh_GIA, loc = rbind(GPS_loc, Mesh_GIA$loc))
  
  ## Create the estimation and prediction stack
  st.est <- inla.stack(data = list(y=GPS_data), A = list(A_data),
                       effects = list(GIA = 1:GIA_spde$n.spde), tag = "est")
  st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                        effects = list(GIA=1:GIA_spde$n.spde), tag = "pred")
  stGIA <- inla.stack(st.est, st.pred)
  
  ## Fix the GPS errors
  hyper <- list(prec = list(fixed = TRUE, initial = 0))
  formula = y ~ -1 +  f(GIA, model = GIA_spde)
  prec_scale <- c(1/GPS$std^2, rep(1, nrow(A_pred)))
  
  ## Run INLA
  res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = GIA_spde), family = "gaussian",
                   scale =prec_scale, control.family = list(hyper = hyper),
                   control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))
  INLA_pred <- res_inla$summary.linear.predictor
  
  ## Extract and project predictions
  pred_idx <- inla.stack.index(stGIA, tag = "pred")$data
  GPS_idx <- pred_idx[1:nrow(GPS)]
  GIA_idx <- pred_idx[-c(1:nrow(GPS))]
  
  ## GPS 
  GPS_u <- INLA_pred$sd[GPS_idx]
  GPS_pred <- data.frame(lon = GPSX, lat = GPS$lat, u = GPS_u)
  
  ## GIA
  GIA_diff <- INLA_pred$mean[GIA_idx] 
  GIA_m <- GIA_diff + GIA_prior
  GIA_u <- INLA_pred$sd[GIA_idx]
  proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(360,180), xlim = c(0, 360), ylim = c(-90, 90))
  GIA_grid <- expand.grid(proj$x, proj$y)
  GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                         diff = as.vector(inla.mesh.project(proj, as.vector(GIA_diff))),
                         mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                         u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))))
  
  return(list(res_inla = res_inla, spde = GIA_spde, st = stGIA, 
              mesh = Mesh_GIA, GPS_pred = GPS_pred, GIA_pred = GIA_pred))
}


## 5 Plot the predicted field
map_res <- function(data, xname, yname, fillvar, colpal = NULL, limits=NULL, title){
  plot_data <- data[c(xname, yname, fillvar)]
  names(plot_data) <- c("X", "Y", "Fill")
  beauty <- 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", colour = 'white'), 
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10), 
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.line = element_line(size = 1),
          plot.title = element_text(hjust = 0.5, size = 15),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          panel.border = element_blank())
  world_map <- map_data("world2")
  baseworld <- geom_polygon(data = world_map, aes(x=long, y=lat, group=group), colour="grey", fill = NA, alpha = 0.5)  
  
  colbar <- guide_colorbar(barwidth = 2, barheight = 10, label.position = "right", title.position = "bottom")
  
  if(is.null(colpal)){
    Map <- ggplot(plot_data) + geom_raster(aes(x = X, y = Y, fill = Fill)) + coord_fixed() + 
      xlab("Longitude") + ylab("Latitude") + 
      scale_x_continuous(limits=c(0,360),  expand = c(0, 0)) + 
      scale_y_continuous(limits=c(-90,90),  expand = c(0, 0)) + 
      scale_fill_gradient2(low = "#0000BF", mid = "white", high = "#BF0000", midpoint = 0, 
                           name = "mm/yr", limits = limits,  guide = colbar) 
  }else{
    Map <- ggplot(plot_data) + geom_raster(aes(x = X, y = Y, fill = Fill)) + coord_fixed() + 
      xlab("Longitude") + ylab("Latitude") + 
      scale_x_continuous(limits=c(0,360),  expand = c(0, 0)) + 
      scale_y_continuous(limits=c(-90,90),  expand = c(0, 0)) + 
      scale_fill_gradientn(colors = colpal, name = "mm/yr", limits = limits,  guide = colbar) 
  }
  
  Map <- Map +  baseworld + ggtitle(title) + beauty
  return(Map)
}


## lattice plot
level_GIA <- function(pred_data, var, style = c("zero", "positive"), layout, title){
  range_var <- range(pred_data[,var], na.rm = TRUE)
  NAid <- which(is.na(pred_data[, var]))
  pred_data[NAid,var] <- range_var[2] + 1
  plotdata <- data.frame(lon = pred_data$lon, lat = pred_data$lat, var = pred_data[, var], model = pred_data$model)
  
  if(style == "zero"){
    brks1 <- seq(floor(range_var[1]), -0.5, length.out=9)
    brks2 <- c(seq(0.5, floor(range_var[2]), length.out = 7), range_var[2]+ 0.9, ceiling(range_var[2] + 1))
    colreg1 <- rev(RColorBrewer::brewer.pal(9, "Blues"))
    colreg2 <- RColorBrewer::brewer.pal(8, "Reds")
    colreg<- c(colreg1, "#FFFFFF", colreg2, "#bebebe")
    brks <- c(brks1, -0.1, 0.1, brks2)
    lab_tik <- seq(-5*(abs(brks[1])%/%5), 5*(brks[length(brks)-1] %/%5), 5)
    labs <- list(labels = c(lab_tik, "NA"), at = c(lab_tik, (brks[length(brks)-1] + brks[length(brks)])*0.5))
  }else if(style == "positive"){
    brks <- c(seq(0, 0.8, 0.1), seq(1,floor(range_var[2])-1, 0.5), floor(range_var[2]), ceiling(range_var[2]), range_var[2]+1)
    colreg <- c(colorRamps::matlab.like(length(brks)-2), "#bebebe")
    lab_tik <- c(0, 1, seq(2, brks[length(brks)-1], 2))
    labs <- list(labels = c(lab_tik, "NA"), at = c(lab_tik, (brks[length(brks)-1] + brks[length(brks)])*0.5))
  }
  
  lattice::levelplot(var ~ lon + lat|model, data = plotdata, aspect = "iso",
                     col.regions= colreg, at=brks, 
                     colorkey = list(col = colreg, labels = labs),
                     panel = function(x,y,z,...){
                       lattice::panel.levelplot(x,y,z,...)
                       map2 <- map("world2", interior = FALSE, plot = FALSE)
                       lattice::panel.xyplot(x=map2$x, y=map2$y, type = "l", col = "black")
                     },
                     layout = layout,
                     main = title, xlab = "longitude", ylab = "latitude")
}



## plot simulated field
map_res2 <- function(res, colpal = NULL, limits=NULL, title, ...){
  mesh <- res$mesh
  proj = inla.mesh.projector(mesh,  projection = "longlat", dims = c(360,180), xlim = c(0, 360), ylim = c(-90, 90))
  field.pred <- inla.mesh.project(proj, as.vector(res$GIA_predn[,field]))
  image.plot(list(x = proj$x, y=proj$y, z = field.pred), ...)
  map("world2", interior = FALSE, add = T, col = "grey")
}



## 6 Plot zoom in
map_zoom <-function(data_field, data_obs, zoom_coords, colpal, limits=NULL){
  lon1 <- zoom_coords$lon[1]
  lon2 <- zoom_coords$lon[2]
  lat1 <- zoom_coords$lat[1]
  lat2 <- zoom_coords$lat[2]
  
  beauty <- 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", colour = 'white'), 
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10), 
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.line = element_line(size = 1),
          plot.title = element_text(hjust = 0.5, size = 15),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          panel.border = element_blank())
  world_map <- map_data("world2")
  baseworld <- geom_polygon(data = world_map, aes(x=long, y=lat, group=group), colour="grey", fill = NA, alpha = 0.5) + beauty 
  
  
  colbar <- guide_colorbar(barwidth = 2, barheight = 10, label.position = "right", title.position = "bottom")
  
  zoom_data <- subset(data_field, lon > lon1 & lon < lon2 & lat > lat1 & lat < lat2)
  zoom_GPS <- subset(data_obs, lon > lon1 & lon < lon2 & lat > lat1 & lat < lat2)
  
  map_zoom <- ggplot(data=zoom_data) + geom_raster(aes(x = lon, y = lat, fill = u)) + 
    coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
    scale_x_continuous(limits=c(lon1,lon2),  expand = c(0.03, 0.03)) + 
    scale_y_continuous(limits=c(lat1,lat2),  expand = c(0, 0)) +
    scale_fill_gradientn(colors = colpal, name = "mm/yr", guide = colbar) 
  
  zoom_std <- subset(zoom_data, lon %in% seq(lon1, lon2, 2))
  zoom_std <- subset(zoom_std, lat %in% seq(lat1, lat2, 2))
  
  map_zoom <- map_zoom + geom_point(data=zoom_std, aes(x=lon, y=lat), pch=19, size = zoom_std$u,
                                    col = "blue", fill = "blue", alpha=0.6) 
  map_zoom <- map_zoom + geom_point(data=zoom_GPS, aes(x=lon, y=lat), pch= "*", size = 10,
                                    col = "red", alpha = 0.8) 
  
  ## Calculate the approximate correlation length
  cor_deg <- round(rho_mode /(2*pi) * 360)
  cor_dist <- round(rho_mode * 6371)
  
  map_zoom <- map_zoom + baseworld  + beauty + 
    ggtitle("GIA uncertainty Zoom-in", 
            subtitle = paste("The estimated correlation length is about", eval(cor_dist), "km", "(", eval(cor_deg),"degree)"))
  
}


## 7 get the subset of the mesh for buiding the Q
mesh.sub <- function(mesh, Omega, i = 2, sphere=TRUE){
  mesh_sub <- mesh
  id_tri <- sort(Omega[[i]]) # triangle id in the subset
  id_vet <- sort(unique(as.vector(mesh$graph$tv[Omega[[i]], ]))) # vertex id in the subset
  
  if(!sphere){
    ## get the subset of the graph
    sub_tv <- apply(mesh$graph$tv[id_tri,], 2, match, table = id_vet)
    sub_vt <- match(mesh$graph$vt[id_vet,], id_tri)
    sub_tt <- apply(mesh$graph$tt[id_tri,], 2, match, table = id_tri)
    sub_tti <- mesh$graph$tti[id_tri,]
    sub_vv <- mesh$graph$vv[id_vet,id_vet]
    
    ## Now assemble 
    mesh_sub$t <- length(id_tri) 
    mesh_sub$n <- length(id_vet)
    mesh_sub$loc <- mesh$loc[id_vet,]
    mesh_sub$graph$tv <- sub_tv
    mesh_sub$graph$vt <- sub_vt
    mesh_sub$graph$tt <- sub_tt
    mesh_sub$graph$tti <- sub_tti
    mesh_sub$graph$vv <- sub_vv
    
    mesh_sub$idx$loc <- na.omit(match(mesh$idx$loc, id_vet))
    mesh_sub$idx$segm <- na.omit(match(mesh$idx$segm, id_vet))
    mesh_sub$segm$bnd$idx <- na.omit(apply(mesh$segm$bnd$idx, 2, match, table = id_vet))
    
    mesh_sub$posTri <- mesh_sub$posTri[Omega[[i]],]
    mesh_sub$Trill <- mesh_sub$Trill[Omega[[i]],]
    mesh_sub$Vidx <- id_vet # the vertice id in the original mesh
  }
  else{
    ## Now update the index for all
    id_all_v <- mesh$idx$loc
    id_all_v[-id_vet] <- NA
    id_all_v[which(!is.na(id_all_v))] <- 1:length(id_vet)
    
    id_all_t <- 1:mesh$t
    id_all_t[-id_tri] <- NA
    id_all_t[which(!is.na(id_all_t))] <- 1:length(id_tri)
    ## Modify the mesh to be the subset
    
    
    ## get the subset of the graph
    sub_tv <- mesh$graph$tv[id_tri, ]
    sub_vt <- mesh$graph$vt[id_vet,]
    sub_tt <- mesh$graph$tt[id_tri,]
    sub_tti <- mesh$graph$tti[id_tri,]
    sub_vv <- mesh$graph$vv[id_vet,id_vet]
    
    ## Change the triangle indices of the subset
    sub_tv2 <- cbind(id_all_v[sub_tv[,1]], id_all_v[sub_tv[,2]], id_all_v[sub_tv[,3]])
    sub_vt2 <- id_all_t[sub_vt]
    sub_tt2 <- cbind(id_all_t[sub_tt[,1]], id_all_t[sub_tt[,2]], id_all_t[sub_tt[,3]])
    
    ## Now assemble 
    mesh_sub$t <- length(id_tri) 
    mesh_sub$n <- length(id_vet)
    mesh_sub$loc <- mesh_sub$loc[sort(id_vet), ]
    mesh_sub$graph$tv <- sub_tv2
    mesh_sub$graph$vt <- sub_vt2
    mesh_sub$graph$tt <- sub_tt2
    mesh_sub$graph$tti <- sub_tti
    mesh_sub$graph$vv <- sub_vv
    
    mesh_sub$idx$loc <- 1:length(id_vet)
    mesh_sub$posTri <- mesh_sub$posTri[Omega[[i]],]
    mesh_sub$Trill <- mesh_sub$Trill[Omega[[i]],]
    mesh_sub$Vidx <- id_vet # the vertice id in the original mesh
  }
  return(mesh_sub)
}

#### additional plot function
## plot posterior hyperparameters
marginal_par <- function(res, process, partition = FALSE, theta.names=NULL, plot = FALSE){
  if(partition){
    if(is.null(theta.names)) stop("theta.names cannot be NULL")
    theta.list <- list()
    ntheta <- length(theta.names)
    theta.modes <- rep(0, ntheta)
    res_inla <- res$res_inla
    for (i in 1:ntheta){
      thetai <- res_inla$marginals.hyperpar[[i]]
      theta.list[[i]]<- inla.tmarginal(exp, thetai)
      ltheta_mode <- res_inla$summary.hyperpar$mode[i]
      ltheta_mean <- res_inla$summary.hyperpar$mean[i]
      ltheta_sd <- res_inla$summary.hyperpar$sd[i]
      theta.modes[i] <- exp(ltheta_mean - ltheta_sd^2)
    }
    names(theta.list) <- theta.names
    if(plot){
      for(i in 1:ntheta)
        plot(theta.list[[i]], type = "l", main = bquote(bold(theta("mode") == .(round(theta.modes[i], 4)))))
    }
    return(list(thetamars = theta.list, thetam = theta.modes))
  }else{
    res_inla <- res$res_inla
    spde <- res$spde[[process]]
    pars <- inla.spde2.result(res_inla, process, spde, do.transf=TRUE)
    Vmar <-pars$marginals.variance.nominal[[1]]
    Rmar <- pars$marginals.range.nominal[[1]]
    
    ## Find the mode of rho and sigma^2
    lrho_mode <- pars$summary.log.range.nominal$mode
    lrho_mean <- pars$summary.log.range.nominal$mean
    lrho_sd <- pars$summary.log.range.nominal$sd
    rho_mode <- exp(lrho_mean - lrho_sd^2)
    
    lsigma_mode <- pars$summary.log.variance.nominal$mode
    lsigma_mean <- pars$summary.log.variance.nominal$mean
    lsigma_sd <- pars$summary.log.variance.nominal$sd
    sigma_mode <- exp(lsigma_mean - lsigma_sd^2)
    if(plot){
      plot(Vmar, type = "l", main = bquote(bold({sigma^2}("mode") == .(round(sigma_mode, 4)))))
      plot(Rmar, type = "l", main = bquote(bold(rho("mode") == .(round(rho_mode, 4)))))
    }
    return(list(Vmar=Vmar, Rmar=Rmar, sigma_mode=sigma_mode, rho_mode=rho_mode))
  }
}






