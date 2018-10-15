## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))

## The formular -- we add the constraint that mass change sum to zero
## constraint 1 -- vertices sum to zero
A1 <- matrix(1, nrow = 1, ncol = mesh0$n)
## Constraint 2 -- grace loc sum to zero
A2 <- matrix(colSums(A_grace), nrow = 1)
## Constraint 3 -- longlat grid sum to zero
gridll <- gia_sp@coords
gridxyz <- do.call(cbind, Lll2xyz(lat = gridll[,2], lon = gridll[,1]))
A_grid <- inla.spde.make.A(mesh = mesh0, loc = gridxyz)
weights <- gia_sp$areas/sum(gia_sp$areas)
A3 <- matrix(weights, nrow = 1) %*% A_grid
A <- as.matrix(rbind(A1, A2, A3))
## The formular -- we add the constraint that mass change sum to zero
formula = y ~ -1 + f(mass, model = mass_spde, extraconstr = list(A = A, e = c(0, 0, 0))) +  f(steric, model = steric_spde)

if(errorname == "raw"){
  prec_scale <- prec_scale
  warning("regional data has no raw error; use inflated error instead.")
}else if(errorname == "inflate"){
  prec_scale <- prec_scale
}else{
  prec_scale <- prec_scale1
}

## Run INLA
## No GIA error
#INLA:::inla.dynload.workaround() 
res_inla <- inla(formula, data = inla.stack.data(stkall), family = "gaussian",
                 scale =prec_scale, 
                 control.family = list(hyper = hyper), control.compute = list(config = TRUE),
                 control.predictor=list(A=inla.stack.A(stkall), compute = TRUE), 
                 verbose = TRUE)
saveRDS(res_inla, file = paste0(dd_root, "WP1-BHM/Experiment2/", filename, "_inla.rds"))
