###### This script run the inla procedure for Experiment 3
####################################################################################################

prec_scale <- c(1/grace_northA$std^2,  1/GPS$std^2)

## Fix observation errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))


## The formular -- we add the constraint that mass change sum to zero
formula = y ~ -1 + f(mass, model = mass_spde) +  f(GIA, model = gia_spde)

res_inla <- inla(formula, data = inla.stack.data(stkall), family = "gaussian",
                 scale =prec_scale, 
                 control.family = list(hyper = hyper), control.compute = list(config = TRUE),
                 control.predictor=list(A=inla.stack.A(stkall), compute = TRUE), 
                 verbose = TRUE)

saveRDS(res_inla, file = paste0(dd_save,expname, "_inla.rds"))
