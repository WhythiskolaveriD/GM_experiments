####################################
#### Experiment 1b Run INLA     ####
#### plug-in approximation      ####
####################################
xx <- commandArgs(TRUE)
eval(parse(text=xx))
source(infile)

#### 0 path and file setup
wkdir <- getwd()
outdir <- "/./projects/GlobalMass/WP1-BHM/Experiment1b/output/"
#.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.3/")
data_path <- "/./projects/GlobalMass/WP1-BHM/Experiment1b/Data_inputs/"

GIA_name <- strsplit(strsplit(GIA_input, split = ".txt")[[1]], split="_")[[1]][2]
GPS_name <- strsplit(strsplit(GPS_input, split = ".txt")[[1]], split="_")[[1]][2]

GIA_file <- paste0(data_path, "GIA/", GIA_input)
GPS_file <- paste0(data_path, "GPS/", GPS_input)

outname <- paste("inla", GIA_name, GPS_name, sep="_" ) 
  
  
#### 1 Prep Priors 
## Priors mean and variance for rho and sigma
mu_r <- 500/6371
v_r <- (1000/6371)^2
mu_s <- 20
v_s <- 40^2

source("~/glbm/Experiment1b/Rscript/1Prep_data.R")
source("~/glbm/Experiment1b/Rscript/2INLA_infer.R")
source("~/glbm/Experiment1b/Rscript/3plot_INLA.R")
