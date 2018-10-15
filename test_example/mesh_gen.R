#### This file contains code for generating different meshes on sphere


#### Generate equally distributed points around a sphere -- Fibonacci spiral
fiboSphere1 <- function(n = 1000L) {
## reference:
  ## 1 http://blog.marmakoide.org/?p=1
  ## 2 Measurement of Areas on a Sphere Using Fibonacci and Latitude–Longitude Lattices (2010)
  ## 3 http://lgdv.cs.fau.de/uploads/publications/spherical_fibonacci_mapping_opt.pdf
  phi <- (sqrt(5) + 1) / 2  # golden ratio
  i <- seq(-(n-1), (n-1),2)
  theta <- 2 * pi * i / phi
  z <- i/n
  radius <- sqrt((1-z^2))
  
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  return(cbind(x,y,z))
}

#### Generate equally distributed points around a sphere -- Fibonacci spiral
fiboSphere2 <- function(N = 1000L, L0 = FALSE) {
  ## Reference (note that points generated from 2D are slightly different from 3D)
  ## Measurement of Areas on a Sphere Using Fibonacci and Latitude–Longitude Lattices (2010)
  phi <- (sqrt(5) + 1) / 2 - 1 # golden ratio
  ga <- phi * 2 * pi           # golden angle
  
  i <- seq(-N, N)
  P <- 2 * N + 1
  lat <- asin(2*i / P) * 180 / pi
  if(L0){
  lon <- ((2 * pi * i / phi) %% pi) * 360 / pi
  }else{
    lon <- ((2 * pi * i / phi) %% pi) * 360 / pi - 180
    }
  cbind(lon = lon, lat = lat)
}
library(GEOmap)

aa <- fiboSphere2(1000)
bb <- do.call(cbind, Lll2xyz(lat = aa[,2], lon = aa[,1]))

library(INLA)
Mesh_b <- inla.mesh.2d(loc = bb, cutoff = 0.001, max.edge = 5) ## small

plot(Mesh_b, rgl = T)
points3d(bb, add = T, col = "green", size = 8, pch = 20)
