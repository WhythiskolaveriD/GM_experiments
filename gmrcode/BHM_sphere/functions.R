## Some functions used in the experiments.

## Transforming priors from original space to log-normal space
Tlognorm <- function(mu, v){
  logv <- log(1 + v/mu^2)
  logmu <- log(mu^2) - 0.5*log(mu^2 + v)
  return(c(logmu, logv))
}

## Generate Fibonacci points on the sphere
fiboSphere <- function(N = 1000L, LL = TRUE, L0 = FALSE) {
  ## Reference (note that points generated from 2D are slightly different from 3D)
  ## Measurement of Areas on a Sphere Using Fibonacci and Latitude–Longitude Lattices (2010)
  phi <- (sqrt(5) + 1) / 2  # golden ratio
  if(LL){ ## for lonlat coords
    i <- seq(-N, N)
    P <- 2 * N + 1
    lat <- asin(2*i / P) * 180 / pi
    if(L0){ ## for lon in [0, 360]
      lon <- ((2 * pi * i / phi) %% pi) * 360 / pi
    }else{ ## for lon in [-180, 180]
      lon <- ((2 * pi * i / phi) %% pi) * 360 / pi - 180
    }
    return(cbind(lon = lon, lat = lat))
  }
  else{ ## for xyz coords
    i <- seq(-(N-1), (N-1),2)
    theta <- 2 * pi * i / phi
    sphi <- i/N
    cphi <- sqrt((N+i) * (N-i))/N
    
    x <- cphi * sin(theta)
    y <- cphi * cos(theta)
    z <- sphi
    return(cbind(x,y,z))
  }
}


## Output dataframe as an netcdf file
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

