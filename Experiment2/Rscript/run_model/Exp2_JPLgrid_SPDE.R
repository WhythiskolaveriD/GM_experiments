## Build the SPDEs for the latent processes

## generate the mesh
#########################################################################################################
cache_mesh <- paste0(dd_root,"WP1-BHM/Experiment2/data/exp2_mesh.rda")

if(file.exists(cache_mesh)){
  load(cache_mesh)
}else{
  ## Genereate Fibonacci points on the sphere
  fibo_points <- fiboSphere(N = 3e4, LL = FALSE)
  mesh0 <- inla.mesh.2d(loc = fibo_points, cutoff = 0.01, max.edge = 1)
  summary(mesh0)
  
  ## Identify vertice zone: land (ice or hydrology), boundary, ocean
  ## sep land and ocean vertices, identify verticies in boundary triangles, asign them to be small values!
  mesh_loc <- mesh0$loc
  mesh_locll <- Lxyz2ll(list(x = mesh_loc[,1], y = mesh_loc[,2], z = mesh_loc[,3]))
  mesh_locll <- cbind(mesh_locll$lon, mesh_locll$lat)
  ocean_in <- over(SpatialPoints(coords = mesh_locll, proj4string = CRS(proj4string(Ocean))), as(Ocean, "SpatialPolygons"))
  ocean_id <- which(ocean_in == 2)
  land_id <- (1:length(ocean_in))[-ocean_id]
  
  ## Identify vertices in triangles across boundaries 
  ocean_tri <- apply(mesh0$graph$tv, 2, function(x) match(x, ocean_id)) 
  ocean_trivn <- rowSums(!is.na(ocean_tri))
  ocean_trid <- which(ocean_trivn == 3)
  
  land_tri <-   apply(mesh0$graph$tv, 2, function(x) match(x, land_id))
  land_trivn <- rowSums(!is.na(land_tri))
  land_trid <- which(land_trivn == 3)
  
  bound_trid <- which(land_trivn < 3 & land_trivn > 0) ## checked same when use ocean
  bound_id <- unique(c(mesh0$graph$tv[bound_trid,]))
  #plot(mesh_locll[bound_vid, ])
  
  ## remove boundary points in land and ocean
  ocean_id <- ocean_id[!(ocean_id %in% bound_id)]
  land_id <- land_id[!(land_id %in% bound_id)]
  region_vid <- list(ocean = ocean_id, land = land_id, bound = bound_id)
  
  ## Create the list for devided regions
  region_tid <- list(ocean_trid, land_trid, bound_trid)
  mesh0$t <- nrow(mesh0$graph$tv)
  mesh_ocean <- mesh.sub(mesh0, region_tid, 1)
  mesh_land <- mesh.sub(mesh0, region_tid, 2)
  mesh_bound <- mesh.sub(mesh0, region_tid, 3)
  
  save(mesh0, mesh_locll, region_vid, region_tid, mesh_land, mesh_ocean, mesh_bound, file = cache_mesh)
}


#### Set up the priors
#########################################################################################################
## Priors for mass
# -- for the land
m_rg_land <- 700/6371 # mean of rho
m_sg_land <- 25 # mean of sigma

# -- for the ocean
m_rg_ocean <- 1500/6371 # mean of rho
m_sg_ocean <- 5 # mean of sigma

# -- for boundary
m_rg_bound <- 200/6371
m_sg_bound <- 40

## Priors for steric
mu_r2 <- 2000/6371 # mean of rho
mu_s2 <- 15 # mean of sigma


## Transform the parameters for a lognormal distribution
#########################################################################################################
## -- mass
lsig_landm <- log(m_sg_land)
lsig_lands <- 1 
lrho_landm <- log(m_rg_land)
lrho_lands <- abs(lrho_landm)/2
lkap_land <- log(8)/2 - lrho_landm
ltau_land <- 0.5*log(1/(4*pi)) - lsig_landm - lkap_land

lsig_oceanm <- log(m_sg_ocean)
lsig_oceans <- 1
lrho_oceanm <- log(m_rg_ocean)
lrho_oceans <- abs(lrho_oceanm)/2
lkap_ocean <- log(8)/2 - lrho_oceanm
ltau_ocean <- 0.5*log(1/(4*pi)) - lsig_oceanm - lkap_ocean

trho_bound <- log(m_rg_bound)
tsig_bound <- log(m_sg_bound)
lkap_bound <- log(8)/2 - trho_bound
ltau_bound <- 0.5*log(1/(4*pi)) - tsig_bound - lkap_bound

## -- steric
lsigma2 <- log(mu_s2)
theta2_ss <- 1
lrho2 <- log(mu_r2)
theta2_rs <- abs(lrho2)/2
lkappa2 <- log(8)/2 - lrho2
ltau2 <- 0.5*log(1/(4*pi)) - lsigma2 - lkappa2


#### Build the SPDEs
#########################################################################################################
## Build the mass spde
## Define the basis functions
## define basis function for tau
nv <- nrow(mesh0$loc)
np <- 5 
Btmat <- Bkmat <- matrix(0, nrow = np, ncol = nv)
## ocean basis
bt1 <- c(ltau_ocean, -1, 1, 0, 0)
bk1 <- c(lkap_ocean, 0, -1, 0, 0)
## land basis
bt2 <- c(ltau_land, 0, 0, -1, 1)
bk2 <- c(lkap_land, 0, 0, 0, -1)
## boundary basis
btb <- c(ltau_bound, 0, 0, 0, 0)
bkb <- c(lkap_bound, 0, 0, 0, 0)

Btmat[,region_vid$ocean] <- bt1
Btmat[,region_vid$land] <- bt2
Btmat[,region_vid$bound] <- btb
Btmat <- t(Btmat)

Bkmat[,region_vid$ocean] <- bk1
Bkmat[,region_vid$land] <- bk2
Bkmat[,region_vid$bound] <- bkb
Bkmat <- t(Bkmat)


mass_spde <- inla.spde2.matern(mesh0, B.tau = Btmat, B.kappa = Bkmat, theta.prior.mean = c(0,0,0,0), 
                               theta.prior.prec = c(sqrt(1/lsig_oceans), sqrt(1/lrho_oceans), 
                                                    sqrt(1/lsig_lands), sqrt(1/lrho_lands)))

## Build the steric spde
steric_spde <- inla.spde2.matern(mesh_ocean, B.tau = matrix(c(ltau2, -1, 1),1,3), B.kappa = matrix(c(lkappa2, 0, -1), 1,3),
                                 theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta2_ss), sqrt(1/theta2_rs)))






