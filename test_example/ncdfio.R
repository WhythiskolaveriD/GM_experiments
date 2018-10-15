#### This file contains functions for i/o ncdf

#### save a dataframe into an ncdf
#### depend on library ncdf4
library(ncdf4)

df2ncdf <- function(df, fname, vars, title = NULL, append = FALSE){
  # - df: a dataframe object that contains the coordinates and variables to be into a ncdf file
  #       It should contain coordinates with name "lon" and "lat".
  #       The dataframe should be ordered by first by lat and then lon
  # - fname: a charater string for the file name
  # - vars: a list of character vectors and each element contains the infomation
  #         -- the variable name in the df
  #         -- the short name of the variable
  #         -- the long name of the variable
  #         -- the unit
  # - append: add new variable values in existing ncdf file
  
  ## define dimensions
  lons <- unique(df$lon)
  lats <- unique(df$lat)
  llon <- length(lons)
  llat <- length(lats)
  londim <- ncdim_def("lon", "degrees_east", as.double(lons))
  latdim <- ncdim_def("lat", "degrees_north", as.double(lats))
  
  ## define variables
  nvar <- ncol(df) - 2
  var_def <- lapply(vars, function(x) ncvar_def(name = x[2], units = x[4], dim = list(londim, latdim), missval = NA, 
                                                longname=x[3], prec = "single"))
  
  ## create a netCDF file 
  if(append){
    nc_open(fname)
    ## put the array into the file
    for(i in 1:nvar){
      ncvar_put(ncout, var_def[i], array(df[,vars[[i]][1]], dim = c(llon, llat)))
    }
  }else{
    ncout <- nc_create(fname, var_def, force_v4 = TRUE)
    ## put the array into the file
    for(i in 1:nvar){
      ncvar_put(ncout, var_def[[i]], array(df[,vars[[i]][1]], dim = c(llon, llat)))
    }
    
    ## put additional attributes into dimension and data variables
    ncatt_put(ncout, "lon", "axis", "X")  
    ncatt_put(ncout, "lat", "axis", "Y")
    
    ## add global attributes
    ncatt_put(ncout, 0, "title", title)
    }
  
  # close the file, writing data to disk
  nc_close(ncout)
}

# ## A test using steri_bhm
# df2ncdf(df = dftest, fname = "test.nc",
#         vars = list(c("trend", "trend", "predicted steric trend", "mm/yr"), 
#                                  c("err", "err", "predicted uncertainty", "mm/yr")),
#         title = "bhm steric")
