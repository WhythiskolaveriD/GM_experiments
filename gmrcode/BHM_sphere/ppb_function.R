####################### DESCRIPTION #######################
###     Functions for the parameter partition model     ###
###                    with buffer                      ###
###     Defined based on inla.doc("rgeneric")           ###
###########################################################


### LOAD LIBRARIES
library(INLA)
library(rgeos)
library(fields)
library(rgdal)

############################################################
### DEFINE THE PROCESS PARTITION MODEL
ppb.inla.model = function(Q, log.prior) {
  ## Input
  # Q and log.prior
  # - must be functions of theta
  # - must be self-contained, with constants fixed in environment
  
  ## Input verification
  if(identical(environment(Q),.GlobalEnv)) stop("Cannot use global environments")
  if(identical(environment(log.prior),.GlobalEnv)) stop("Cannot use global environments")
  if(anyDuplicated(ls(environment(Q)), ls(environment(log.prior)))) stop("Duplicate variables")
  
  ## Define th emodel
  model.rgeneric = inla.rgeneric.define(model = ppb.rgeneric.model, Q = Q, log.prior = log.prior)
  return(model.rgeneric)
}

############################################################
### DEFINE THE MODEL SKELETON FUNCTION
ppb.rgeneric.model = function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), 
  theta = NULL)
{
  
  ## Define the Graph matrix (non zero elements in Q)
  graph = function(theta)
  {
    ## Note that Q is defined by its lower triangle elements
    require(methods)
    ntheta = length(initial())
    
    G1 = Q(theta=(1:ntheta)/3.217233456)
    G1[G1 != 0] = 1
    G2 = Q(theta=(1:ntheta)^2/12.1543534)
    G2[G2 != 0] = 1
    
    return (G1+G2)
    # - this should never give zeroes by random.
  }
  
  ## Q is given outside at the top level function
  
  ## Define the mean of the model
  mu <- function(theta) numeric(0)
  
  ## Define the initial values for theta
  initial = function(theta)
  {
    return (environment(Q)$initial.theta)
  }
  
  ## Define the noarmalising consant -- let INLA compute it from Q
  log.norm.const = function(theta) numeric(0)
  
  ## The log.prior is given outside at the top level function
  
  ## exiting the C-program
  quit <- function(theta) invisible()
  
  ## Retrun values
  val = do.call(match.arg(cmd), args = list(theta))
  return (val)
}
############################################################
### BUILD THE FINITE ELEMENT MATRICES
ppb.fem <- function(mesh, subdomain, sphere){
  ## input:
  # - mesh: the entire mesh
  # - subdomain: list of triangle indices in each subdomain
  
  p <- length(subdomain)
  
  #fem_all <- inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2, output = list("c0", "c1", "g1", "g2"))
  vn <- mesh$n
  femList <- list()
  for (i in 1:p){
    Cmat <- G1 <- G2 <- Matrix(0, vn, vn)
    meshsub <- mesh.sub(mesh, subdomain, i, sphere)
    femsub <- inla.fmesher.smorg(meshsub$loc, meshsub$graph$tv, fem =2, output = list("c0", "c1", "g1", "g2"))
    node_id <- meshsub$Vidx
    
    Cmat[node_id, node_id] <- femsub$c0
    G1[node_id, node_id] <- femsub$g1
    G2[node_id, node_id] <- femsub$g2
    Cmat <- tril(Cmat)
    G1 <- tril(G1)
    G2 <- tril(G2)
    
    femi <- list(Cmat = Cmat, G1 = G1, G2 = G2)
    femList[[i]] <- femi
  }
  return(femList)
}


############################################################
### BUILD THE PRECISION MATRIX Q
ppb.precision <- function(femlist, theta, theta_buffer){
  ## Assume an spde2 model, so alpha = 2
  ## The model is on one R2 and S2, so nv = 1
  
  ## Input:
  # - fems: a list of all the finite elment matrics (lower triangular) for each subset mesh
  #         the fems of buffer zone should come in the end!
  # - theta:  a (p-1) x 2 lenght vector for the parameters, p is the number of independent subsets
  #           the first half are log(rho) and the second are log(sigma)
  # - theta_buffer: a vector for log(rho) and log(sigma) in the buffer zone, must be given
  
  ## check if theta_buffer is null
  if(is.null(theta_buffer)) stop("theta_buffer must not be NULL!")
  
  p <- length(femlist)
  lrhos <- c(theta[1:(p-1)], theta_buffer[1])
  lsigmas <- c(theta[-(1:(p-1))], theta_buffer[2])
  
  
  ## For a subset i, use the INLA spde functions to build the Q_i
  vn <- nrow(femlist[[1]]$Cmat)
  Q <- Matrix(0, vn, vn)
  Tau  <- rep(0, vn)
  
  for (i in 1:p){
    Kappa <- exp(0.5*log(8) - lrhos[i])
    Tau[femlist[[i]]$node_id] <- exp(-0.5*log(4*pi) - lsigmas[i] - log(Kappa))
    femi <- femlist[[i]]
    Q <- Q + (Kappa^4*femi$Cmat + 2*Kappa^2*femi$G1 +femi$G2)
  }
  Tau <- Diagonal(n = vn, x = Tau)
  Q <- Tau %*% Q %*% Tau
  ## Assemble these Qi to Q
  Q <- inla.as.dgTMatrix(Q)
  return(Q)
  return(Q)
}



############################################################
## DEFINE Q FUNCTION IN THE SKELETON
ppb.create.Q = function(mesh, subdomain, sphere, initial.theta = NULL, theta_buffer)
{
  ## Input:
  # - mesh: the entire mesh
  # - initial.theta: initial values for theta, can be null
  # - subdomain: list of triangle indices in each subdomain
  ppb.precision = ppb.precision
  p <- length(subdomain)
  if(is.null(initial.theta)) {initial.theta = rep(0, p*2)}
  
  femL <- ppb.fem(mesh = mesh, subdomain = subdomain, sphere = sphere)
  theta_buffer <- theta_buffer
  ## Create Q function  
  Q <- function(theta){
    return(ppb.precision(femlist = femL, theta = theta, theta_buffer = theta_buffer))
  }  
  return(Q)
}




############################################################
### DEFINE THE LOG.PRIOR FUNCTION
ppb.create.prior.log.norm = function(prior.param) {
  ## Create the log normal prior for sigma and range in theta
  ## range and sigma follows log normal distribution 
  ## prior.param are the mean and variance of theta
  ## prior.param$sigma[1,1] = E(sigma_1), prior.param$sigma[1,2] = V(sigma_1), ...
  ## prior.param$rho[1,1] = E(rho_1), prior.param$rho[1,2] = V(rho_1), ...
  ## theta is a vector of p*2 with first p elements log(rho) and second p log(sigma)
  
  ## Move to current scope (environment)
  prior.param = prior.param
  
  log.prior = function(theta) {
    rhos_m <- as.numeric(prior.param$rho[,1])
    rhos_s <- as.numeric(prior.param$rho[,2])
    sigmas_m <- as.numeric(prior.param$sigma[,1])
    sigmas_s <- as.numeric(prior.param$sigma[,2])
    
    p <- length(theta)/2
    lrhos <- theta[1:p]
    lsigmas <- theta[-c(1:p)]
    
    val <- 0
    ## Prior for rhos
    for (i in 1:p){
      val <- val + dlnorm(exp(lrhos[i]), meanlog = rhos_m[i], sd = rhos_s[i], log = TRUE)
    }
    ## Prior for sigams
    for (i in 1:p) {
      val <- val + dlnorm(exp(lsigmas[i]), meanlog = sigmas_m[i], sd = sigmas_s[i], log = TRUE)
    }
    return(val)
  }
  return(log.prior)
  # - this environment includes the prior parameters
}

############################################################
### DEFINE THE LOG.PRIOR FUNCTION
ppb.create.priors2 = function(prior.param) {
  ## Create the beta distribution for the range and gamma for sigma
  ## prior.param$sigma[1,1] = alpha, prior.param$sigma[1,2] = beta in the gamma distribution, ...
  ## prior.param$rho[1,1] = alpha, prior.param$rho[1,2] = beta in the beta distribution, ...
  ## theta is a vector of p*2 with first p elements log(rho) and second p log(sigma)
  
  ## Move to current scope (environment)
  prior.param = prior.param
  
  log.prior = function(theta) {
    rhos_a <- as.numeric(prior.param$rho[,1])
    rhos_b <- as.numeric(prior.param$rho[,2])
    sigmas_a <- as.numeric(prior.param$sigma[,1])
    sigmas_b <- as.numeric(prior.param$sigma[,2])
    
    p <- length(theta)/2
    lrhos <- theta[1:p]
    lsigmas <- theta[-c(1:p)]
    
    val <- 0
    ## Prior for rhos
    for (i in 1:p){
      val <- val + dbeta(exp(lrhos[i]),  shape1 = rhos_a[i], shape2 = rhos_b[i], log = TRUE)
    }
    ## Prior for sigams
    for (i in 1:p) {
      val <- val + dgamma(exp(lsigmas[i]), shape = sigmas_a[i], rate = sigmas_b[i], log = TRUE)
    }
    return(val)
  }
  return(log.prior)
  # - this environment includes the prior parameters
}


############################################################
### PLOT THE POSTERIOR PDF OF PARAMETERS
marginal_par <- function(res, process, SPDE2 = FALSE, theta.names=NULL, plot = FALSE){
  if(!SPDE2){
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
        plot(theta.list[[i]], type = "l", main = bquote(bold(.(theta.names[i])("mode") == .(round(theta.modes[i], 4)))))
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
