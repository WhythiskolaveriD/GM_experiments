###### This script build the SPDE and stacks for Experiment 3
####################################################################################################

## Priors for the hyperparameters
## -- mass
mass_rho <- 700/6371
mass_sigma <- 20
lsig_massm <- log(mass_sigma)
lrho_massm <- log(mass_rho)
lkap_mass <- log(8)/2 - lrho_massm
ltau_mass <- 0.5*log(1/(4*pi)) - lsig_massm - lkap_mass
mass_s1 <- 1 
mass_s2 <- abs(lrho_massm)/2

mass_spde <- inla.spde2.matern(mesh0, B.tau = matrix(c(ltau_mass, -1, 1),1,3), B.kappa = matrix(c(lkap_mass, 0, -1), 1,3),
                               theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/mass_s2), sqrt(1/mass_s1)))

## -- GIA: with ice6g as prior mean
## The following gia process is the discrepency process between gia and 
gia_rho <- 1000/6371
gia_sigma <- 10
lsig_giam <- log(gia_sigma)
lrho_giam <- log(gia_rho)
lkap_gia <- log(8)/2 - lrho_giam
ltau_gia <- 0.5*log(1/(4*pi)) - lsig_giam - lkap_gia
gia_s1 <- 1
gia_s2 <- abs(lrho_giam)/2

gia_spde <- inla.spde2.matern(mesh0, B.tau = matrix(c(ltau_gia, -1, 1),1,3), B.kappa = matrix(c(lkap_gia, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/gia_s2), sqrt(1/gia_s1)))


## Link GRACE
ygrace <- grace_northA$detrend

grace_loc <- do.call(cbind, GEOmap::Lll2xyz(lat = grace_northA@coords[,2], lon = grace_northA@coords[,1]))

A_grace2Mass <- inla.spde.make.A(mesh = mesh0, loc = grace_loc)
A_grace2Gia <- inla.spde.make.A(mesh = mesh0, loc = grace_loc)
stk_grace <- inla.stack(data = list(y=ygrace), A = list(A_grace2Mass,  A_grace2Gia),
                        effects = list(list(mass = 1:mass_spde$n.spde),
                                       list(GIA = 1:gia_spde$n.spde)), tag = "grace")

## Link GPS
yGPS <-GPS$detrend
GPS_loc <- do.call(cbind, GEOmap::Lll2xyz(lat = GPS@coords[,2], lon = GPS@coords[,1]))
A_GPS2Gia <- inla.spde.make.A(mesh = mesh0, loc = GPS_loc)
stk_GPS <- inla.stack(data = list(y=yGPS), A = list(A_GPS2Gia),
                      effects = list(list(GIA = 1:gia_spde$n.spde)), tag = "GIA")

stkall <- inla.stack(stk_grace, stk_GPS)