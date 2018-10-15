#### Generate mesh one the sphere except in one polygon
#### Use the world coase line
library(INLA)
library(rgdal)
library(maptools)
library(GEOmap)
library(rgl)

### 1 INLA mesh -- boundary and interior
loc.bnd = matrix(c(0,0, 1,0, 1,1, 0,1), 4, 2, byrow=TRUE)
loc.int = matrix(c(0.3,0.3, 0.8,0.4, 0.7,0.9, 0.2, 0.8, 0.3,0.3), 5, 2, byrow=TRUE)
segm.bnd = inla.mesh.segment(loc.bnd)
segm.int = inla.mesh.segment(loc.int, is.bnd=FALSE)
int.Poly <- SpatialPolygons(list(Polygons(list(Polygon(loc.int)), ID = "1")))


## Points to be meshed
# Want dense grid outside and sparse inside
# first create uniformaly dense points
loc0 <- as.matrix(expand.grid(seq(0, 1, 0.02), seq(0, 1, 0.02)))
# remove points inside the polygon
loc0_id <- unlist(over(int.Poly, SpatialPoints(loc0), returnList=T))
loc0a <- loc0[-loc0_id,]
# sparse grid for the inside
loc1 <- as.matrix(expand.grid(seq(0, 1, 0.1), seq(0, 1, 0.1)))
loc1_id <- unlist(over(int.Poly, SpatialPoints(loc1), returnList=T))
loc1a <- loc1[loc1_id,]
# combine the two
loc <- rbind(loc0a, loc1a)
#loc = matrix(runif(10*2),10,2)*0.9+0.05
mesh = inla.mesh.2d(loc = loc, boundary = segm.bnd, interior = segm.int, max.edge = 0.5, min.angle =40)
# Note that when using interior here offset and max.edge can only be scalars.
# Bondary defines the boundary of the domain of interest
# interior defines a set of segments inside the domain of interest that are desired to the triagle edges.
plot(mesh)
## From this we can build functions to separate polygons and the holes in it.

### 2 meshes for separate domains
## Use the above example we try to build continuous mesh with different shapes within the two domains
## Learn from the different terrain/ barrier model examples
## Idea: Find the triangles inside the polygons and mask them from the mesh when building up the spde model

## Find which triangles are inside the polygon -- centroid incide the polygon

# a) Find all the centroid of all triangles from the mesh
mesh_center <- function(mesh){ # this is exactly the same as Bakka's function "DT.mesh.addon.posTri()"
  mesh$nt <- nrow(mesh$graph$tv)
  mesh$t_center <- matrix(0, nrow = mesh$nt, ncol = 3)
  for (i in 1:mesh$nt){
    t_coord <- mesh$loc[mesh$graph$tv[i,],]
    mesh$t_center[i,] <- colMeans(t_coord)
  }
  return(mesh)
}

# b) Find triangles inside the polygons
mesh2 <- mesh_center(mesh)
insidePoly <- over(int.Poly, SpatialPoints(mesh2$t_center[,1:2]), returnList=T)
plot(mesh)
points(mesh2$t_center[unlist(insidePoly),1:2], col = "blue")

## Mask 1 
## Do inference as before but mask at prediction -- change vertices inside polygons to be zeros in the output


## Mask 2
## Build different a precision matrix based on two different processes 
## extract the Basis matrix and do the operations using the corresponding kappas and taus 
## build the mixture model with the precision matrix and priors

