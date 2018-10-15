#### Rscript to run exp2c_a on a server
#################################################################################
#### Note: gdal need to be loaded on the server!!!!
#### On atlantis: module load /opt/modulefiles/gdal-2.2.3-x86_64
## Set up the path
wd <- "~/"
dd <- "/./projects/GlobalMass/"
gracedata <- "masscons"

## Load libraries and functions
library(rgdal); library(sp);library(GEOmap)
library(INLA)
library(ncdf4)
source(paste0(wd,"gmrcode/BHM_sphere/functions.R"))
source(paste0(wd, "gmrcode/BHM_sphere/partition_fun.R"))

## Load data and prior
source(paste0(wd, "experiments/Experiment2/Doc/c/Rscripts/exp2c_loadData.R"))

## Load spde and stacks
source(paste0(wd, "experiments/Experiment2/Doc/c/Rscripts/exp2c_a_INLAstacks.R"))

## Run INLA
## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))

## The formular -- we add the constraint that mass change sum to zero
formula = y ~ -1 + f(massO, model = massOcean_spde) + f(massL, model = massLand_spde) + 
  f(steric, model = steric_spde)

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stkall), family = "gaussian",
                 scale =prec_scale, 
                 control.family = list(hyper = hyper), control.compute = list(config = TRUE),
                 control.predictor=list(A=inla.stack.A(stkall), compute =TRUE), 
                 verbose = TRUE)


saveRDS(res_inla, file = "/./projects/GlobalMass/WP1-BHM/Experiment2c/Exp2c_a.rds")

## Assemble and save the results
source(paste0(wd, "experiments/Experiment2/Doc/c/Rscripts/exp2c_a_res.R"))