## ----setup0, include=FALSE-----------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## ----setup---------------------------------------------------------------
wd_root <- "~/GM_experiments/" # this is the experiments folder path
wd_script <- "~/GM_experiments/Experiment3/Rscript/" # this is where the Rscripts are
dd_root <- "/./projects/GlobalMass/" # this is where you have the rdsf path
dd_save <- "/./projects/GlobalMass/WP1-BHM/Experiment3/" # This is where you save the results
 
library(rgdal); library(sp);library(GEOmap)
library(INLA)
library(ncdf4)
source(paste0(wd_root,"gmrcode/BHM_sphere/functions.R"))
source(paste0(wd_root, "gmrcode/BHM_sphere/partition_fun.R"))

expname <- filename <- "Exp3"

## ----data----------------------------------------------------------------
source(paste0(wd_script, "preprocessing/", expname, "_data.R"))

## ----mesh----------------------------------------------------------------
source(paste0(wd_script, "preprocessing/", expname, "_meshdata.R"))

## ----SPDEstks------------------------------------------------------------
source(paste0(wd_script, "run_model/", expname, "_SPDEstks.R"))

## ----runinla-------------------------------------------------------------
source(paste0(wd_script, "run_model/", expname, "_runinla.R"))

## ----results-------------------------------------------------------------
source(paste0(wd_script, "analyze_result/", expname, "_res.R"))

