## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## ----expname-------------------------------------------------------------
wd_root <- "~/GM_experiments/" # this is the experiments folder path
wd_script <- "~/GM_experiments/Experiment2/Rscript/" # this is where the Rscripts are
dd_root <- "/./projects/GlobalMass/" # this is where you have the rdsf path
dd_save <- "/./projects/GlobalMass/WP1-BHM/Experiment2/" # This is where you save the results
 
library(rgdal); library(sp);library(GEOmap)
library(INLA)
library(ncdf4)
source(paste0(wd_root,"gmrcode/BHM_sphere/functions.R"))
source(paste0(wd_root, "gmrcode/BHM_sphere/partition_fun.R"))

expname <- "Exp2_GSFCgrid"
errorname <- "inflate"
filename <- paste0(expname, "_", errorname)

## ----S1_data-------------------------------------------------------------
source(paste0(wd_script, "preprocessing/", expname, "_data.R"))

## ----S1_SPDE-------------------------------------------------------------
source(paste0(wd_script, "run_model/", expname, "_SPDE.R"))

## ----S1_stack------------------------------------------------------------
source(paste0(wd_script, "run_model/", expname, "_stacks.R"))

## ----S1_inlarun----------------------------------------------------------
source(paste0(wd_script, "run_model/", expname, "_runINLA.R"))

## ----S1_result-----------------------------------------------------------
source(paste0(wd_script, "analyze_result/", expname,"_res.R"))

## ----getcodes------------------------------------------------------------
knitr::purl("path2rmd/myrmd.Rmd")

