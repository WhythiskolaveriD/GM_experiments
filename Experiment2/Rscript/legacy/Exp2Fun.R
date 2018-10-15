## 1 Transforming priors from original space to log-normal space
Tlognorm <- function(mu, v){
  logv <- log(1 + v/mu^2)
  logmu <- log(mu^2) - 0.5*log(mu^2 + v)
  return(c(logmu, logv))
}


## Wrapper for INLA 
BayesDA_SSH <- function(SSH, GIA, GRACE, trho, tsigma){
  
  ## Transform the GIA data
  Mesh_GIA <- inla.mesh.create(globe= 64) 
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
  proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(720,362), xlim = c(0, 359), ylim = c(-90, 90))
  GIA_grid <- expand.grid(proj$x, proj$y)
  GIA_pred <- data.frame(lon = GIA_grid[,1], lat = GIA_grid[,2],
                         diff = as.vector(inla.mesh.project(proj, as.vector(GIA_diff))),
                         mean = as.vector(inla.mesh.project(proj, as.vector(GIA_m))),
                         u = as.vector(inla.mesh.project(proj, as.vector(GIA_u))))
  
  return(list(res_inla = res_inla, spde = GIA_spde, st = stGIA, 
              mesh = Mesh_GIA, GPS_pred = GPS_pred, GIA_pred = GIA_pred))
}
