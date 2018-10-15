###################################################################
##          A Bayesian Gaussian Update                           ##
## This script performs a full Bayesian Gaussian Update          ##
## of the GIA processes                                          ##
## Given                                                         ##
## 1. prior assumptions of the parameters                        ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model priors                                   ##
###################################################################
## Pseudo Code of the "Slice within Gibbs"                       ##
##  1. Model                                                     ##
## Y | X, e ~ N(Ay %*% X, I * e)                                 ##
## X - Xm ~ GP (0, Sigma(rho, sigma))                            ##
## Xm = X_GIA -- offset                                          ##
## rho ~ U(a, b)                                                 ##
## sigma, e ~ IG (alpha, beta)                                   ##
##  2. Targets: X|Y, rho|Y, Sigma | y, e|Y                       ##
##  3.Algorithm (no coloring!)                                   ##
## r0, s0, e0 intial values, nsamples                            ##
## At Step i                                                     ##
## (1) Update X(i)|sigma(i-1), rho(i-1), e(i-1), y               ##
## (2) Update rho(i)|X(i),y                                      ##
## (3) Update sigma(i)| X(i), y                                  ##
## (4) Update e(i) | X(i), y                                     ##    
###################################################################
#### Setup

### load packages
library(GEOmap)
library(INLA)
if(runserver){
  INLA:::inla.dynload.workaround()
}#Use this on server with old gcc library
library(slice)
library(actuar)
library(spam)
library(parallel)
source("glbm/Experiment1a/Rscript/MVSTplus.R")

### set constants
GPS_loc <- do.call(cbind, Lll2xyz(lat = GPS_obs$lat, lon = GPS_obs$lon))
GPSX <- ifelse(GPS_obs$lon > 180, GPS_obs$lon-360, GPS_obs$lon)
GPSY <- GPS_obs$lat

Mesh_GIA <- inla.mesh.create(globe = 10)
MlocLL <- Lxyz2ll(list(x=Mesh_GIA$loc[,1], y = Mesh_GIA$loc[,2], z = Mesh_GIA$loc[,3]))
MlocLL$lon <- ifelse(MlocLL$lon < 0, MlocLL$lon + 360, MlocLL$lon)
MlocLL$lon <- ifelse(MlocLL$lon > 359.5, MlocLL$lon - 360, MlocLL$lon)
M_sp <- SpatialPoints(data.frame(lon = MlocLL$lon, lat = MlocLL$lat), proj4string = CRS("+proj=longlat")) #This convert GIA_ice6g a SpatialPointDataFrame

Midx <- over(M_sp, Plist)
GIA_mu <- GIA_ice6g_sp$trend[Midx]
## GIA
GIA_ice6g <- read.table("experimentBHM/GIA_truth.txt", header = T)
GIA_ice6g$x_center <- ifelse(GIA_ice6g$x_center > 180, GIA_ice6g$x_center-360, GIA_ice6g$x_center)
GIA_loc6 <- inla.mesh.map(cbind(GIA_ice6g$x_center, GIA_ice6g$y_center), projection = "longlat")

CMat <- inla.spde.make.A(mesh = Mesh_GIA, loc = GPS_loc)
CMat2 <- inla.spde.make.A(mesh = Mesh_GIA, loc = GIA_loc6)
nObs <- nrow(CMat)
nMesh <- ncol(CMat)
x_mu <- matrix(GIA_mu, nrow = nMesh, ncol = 1)
igshape_new <- igshape0 + nObs/2

### Generate synthetic GPS data
y0 <- CMat %*% x_mu
ydata <- GPS_obs$trend
yerr <- GPS_obs$std
Q_obs <- Matrix(0, nrow = nObs, ncol = nObs, sparse = TRUE)
diag(Q_obs) <- (1/yerr)^2

### Transform hyper parameters 
## GMRF Matern covariance parameters
## sigma ~ lognormal(mu_s, v_s)
## rho ~ lognormal(mu_r, v_r)
## Transform to log scale
## log(sigma) = log(sigma0) + theta_1 ~ normal(lmu_s, lv_s)
## log(rho) = log(rho0) + theta_2 ~ normal(lmu_r, lv_r)
tsigma <- Tlognorm(mu_s, v_s)
trho <- Tlognorm(mu_r, v_r)
lmu_s <- tsigma[1]
lv_s <- tsigma[2] # var(theta_1)
lmu_r <- trho[1]
lv_r <- trho[2] # var(theta_2)
## Priors for the SPDE parameters log(kappa) and log(tau):
## log(kappa) = log(kappa0) - theta_2 ~ normal(log(kappa0), lv_r)
## log(tau) = log(tau0) - theta1 + theta2 ~ norm(log(tau0), lv_s + lv_r)
## couple (rho, sigma) and (kappa, tau) through SPDE connection 
lkappa0 <- log(8)/2 - lmu_r
ltau0 <- 0.5*log(1/(4*pi)) - lmu_s - lkappa0
## Build the SPDE
GIA_spde <- inla.spde2.matern(Mesh_GIA, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/lv_s), sqrt(1/lv_r)))



mcmcGIA <- function(initial_vals){
  ### Initial values
  
  ## Latent GMRF
  theta0 <- initial_vals
  Q_GIA0 <- inla.spde.precision(GIA_spde, theta=theta0)
  
  ### Allocate storage
  x_samp <- matrix(0,numsamples, nMesh)
  theta12_samp <- matrix(0, numsamples,2)
  
  
  ### MCMC alogrithms
  if(sampler == "MH_RW"){
    rwsd <- 1
  }
  
  if(sampler == "slice1"){
    slice_theta1 <- slice_theta2 <- slice(d = 1)
    nlearn <- 100
    lscale1 <- lscale2 <- rep(0, nlearn)
  }
  
  if(sampler == "slice2"){
    slice_theta <- slice(d = 2)
    nlearn <- 100 # learning step for the slice sampler
    lscale <- rep(0, nlearn)
  }
  
  ### Start MCMC simulation
  Q_GIA <- Q_GIA0
  theta_t <- theta0
  mmchol <- summary(chol(as.spam.dgCMatrix(Q_GIA)))
  m <- mm <- 0
  t1 <- proc.time()
  pb <- txtProgressBar(min = 0, max = numsamples, style = 3)
  while (m  < numsamples){
    ### 1 Update the latent process
    Q_new <- as.spam.dgCMatrix(crossprod(CMat, Q_obs) %*% CMat + Q_GIA)
    bt <- crossprod(CMat, Q_obs)%*%ydata + Q_GIA %*% x_mu
    cholQ_new <- chol(Q_new,memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices))
    zm_new <- backsolve(cholQ_new, cbind(rnorm(nMesh), forwardsolve(cholQ_new, bt)))
    x_new <- rowSums(zm_new)
    
  
    ### 2 Update the SPDE parameters
    z_GIA <- x_new - x_mu
    
    if(sampler == "MH_RW"){
      log_postcond_theta <- function(theta){
        Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=theta))
        ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
        as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + crossprod(theta/c(sqrt(lv_s), sqrt(lv_r)))) + ldetQ)
      }
      
      r <- 1
      a <- 0
      while(r > a){
        theta_p <- theta_t + rnorm(2, sd = 0.07)
        a <- exp(log_postcond_theta(theta_p) - log_postcond_theta(theta_t))
        r <- runif(1)
      }
      theta_t <- theta_p
    }
    
    ## bivariate slice sampler
    if(sampler == "slice2"){
      log_postcond_theta <- function(theta){
        Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=theta))
        ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
        as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + crossprod(theta/c(sqrt(lv_s), sqrt(lv_r)))) + ldetQ)
      }
      if(m <= nlearn){
        theta_p <- slice_theta(theta_t, log_postcond_theta, learn = TRUE)
        lscale[m] <- getscales(slice_theta)
      }else{
        theta_p <- slice_theta(theta_t, log_postcond_theta, learn = FALSE)
      }
      theta_t <- theta_p
    }
    
    ## univariate slice sampler
    if(sampler == "slice1"){
      log_postcond_theta1 <- function(theta1){
        Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=c(theta1, theta_t[2])))
        ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
        as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + theta1^2/lv_s) + ldetQ)
      }
      if(m <= nlearn){
        theta1_p <- slice_theta1(theta_t[1], log_postcond_theta1, learn = TRUE)
        lscale1[m] <- getscales(slice_theta1)
      }else{
        theta1_p <- slice_theta1(theta_t[1], log_postcond_theta1, learn = FALSE)
      }
      
      
      log_postcond_theta2 <- function(theta2){
        Q_G <- as.spam.dgCMatrix(inla.spde.precision(GIA_spde, theta=c(theta1_p, theta2)))
        ldetQ <- sum(log(diag(chol(Q_G, memory = list(nnzR= mmchol$nnzR ,nnzcolindices = mmchol$nnzcolindices)))))
        as.numeric(-0.5*(crossprod(z_GIA, Q_G) %*% z_GIA + theta2^2/lv_r) + ldetQ)
      }
      if(m <= nlearn){
        theta2_p <- slice_theta2(theta_t[2], log_postcond_theta2, learn = TRUE)
        lscale2[m] <- getscales(slice_theta2)
      }else{
        theta2_p <- slice_theta2(theta_t[2], log_postcond_theta2, learn = FALSE)
      }
      theta_t <- c(theta1_p, theta2_p)
    }
    
    ### 4 Store samples and new values
    if(mm >= burnin){
      if((mm-burnin)%%thinning == 0){
        m <- m +1
        x_samp[m,] <- x_new
        theta12_samp[m,] <- theta_t
      }
    }
    mm <- mm +1
    Q_GIA <- inla.spde.precision(GIA_spde, theta=theta_t)
    
    ### Print the Progression Bar
    setTxtProgressBar(pb, m)
    
  }
  
  ### close the progress bar on finishing
  close(pb)
  
  t2 <- proc.time()
  ttime <- t2-t1
  
  return(list(time = ttime, samples = list(theta1 = theta12_samp[,1],
                                           theta2 = theta12_samp[,2], xGIA = x_samp)))
  
}

source("glbm/Experiment1a/Rscript/mc_hack.R")

res <- mclapply.hack(ini_vals, mcmcGIA)
save(res, file = paste0(wkdir, exname, ".RData"))
