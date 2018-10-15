################################
####  Test MCMC samplers    ####
####                        ####
################################
#### Load mesh and GPS data
wkdir <- getwd()
.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.3/")
runserver = TRUE
max_edge = 0.25
GPS_obs <- read.table("experimentBHM/GPS_combined_20170530_final.txt", header = T)
load("experimentBHM/GIA_sp_info.RData")
#### MCMC parameters
#set.seed(7)#short chain
set.seed(20)#long chain
numsamples = 1000
burnin = 5000
thinning = 50
sampler = "slice1"
n.chains = 3

#### Set priors for the hyper parameters
### Measurement error sigma_e
## sigma_e ~ InverseGamma(shape = 1, rate = 5e-5)
igshape0 <-  1
igscale0 <- 1

### SPDE parameters rho and sigma
## sigma ~ lognormal(mu_s, v_s)
## rho ~ lognormal(mu_r, v_r)
mu_s <- 20
v_s <- 40^2
mu_r <- 500/6371
v_r <- (1000/6371)^2

exname <- paste0("/experimentBHM/1b_", sampler, "_", as.character(numsamples))
ini_vals <- list(vals1 = c(0, 0), vals2 = c(0.5, -0.5), vals3 = c(-0.5, 0.5))
source("glbm/Experiment1b/Rscript/2GP_Bayes_update_mcmc.R")
source("glbm/Experiment1b/Rscript/3GP_Bayes_plot_mcmc.R")
