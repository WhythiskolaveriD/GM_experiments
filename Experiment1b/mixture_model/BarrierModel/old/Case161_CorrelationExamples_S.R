
### DESCRIPTION ###
# This code gives eaxmples of the Barrier model using R-INLA
# See: www.r-inla.org/barrier-model
# In this code we only look at the (a priori) behaviour of the spatial model, the correlation plots, and how to construct meshes
# The simulation in case 162 uses the mesh1 found below
# This code show examples of increasing difficulty in mesh construction

### MASTER LOCATION ###
# Ignore this! Internal use only!
# This document is either a master or a slave
# Master:
# - dtrcode repository
# Slaves: 
# - Gdrive/public/barriermodel1/Case161_CorrelationExamples_S.R
# Synchronise:
# - phdbakka/TutorialDT/TutorialDT.Rnw (parts only)


### LOADING LIBRARIES AND FUNCTIONS ###
rm(list=ls())
source('functMinimalDT.R')

### EXAMPLE 1: NARROW PASSAGE ###
# This example can be thought of as a bridge, or as a barrier with an opening/gap
# In this example:
# - We define a polygon
# - We use interior boundaries in the inla.mesh.segment to make the boundary consistent for different meshes (different max.edge.length)
# - We show how to verify the result

## User input
max.edge.length = 0.4
# - tested: 0.4
smalldist = 0.2

## Mesh setup
loc1 = matrix(c(0,0, 10,0, 0,10, 10,10), 4, 2, byrow = T)
# - this defines the extent of the mesh
# - these are NOT your data locations
# - in an application, if you want the mesh to depend on your data locations, you may use them here

width = .5
p2 = bakka.square.polygon(xlim=c(-1, 5-smalldist/2), ylim=5+width*c(-.5, .5), ret.SP = T) 
# - left square
p3 = bakka.square.polygon(xlim=c(5+smalldist/2, 11), ylim=5+width*c(-.5, .5), ret.SP = T) 
# - right square
p = SpatialPolygons(c(p2@polygons, p3@polygons))
# - the total barrier area polygon
seg = inla.sp2segment(p)
mesh1 = inla.mesh.2d(loc=loc1, interior = seg, max.e = max.edge.length, offset=1)
# - an inla mesh where the triangle edges follow the barrier area polygon

mesh1 = DT.mesh.addon.posTri(mesh1)
# - adding on the positions of the mesh triangles to the mesh object

barrier1 = over(p, SpatialPoints(mesh1$posTri), returnList=T)
# - checking which mesh triangles are inside the barrier area
# - here you may want to use additional polygons (if your later verification shows that you do not have the appropriate model)
barrier1 = unlist(barrier1)
Omega1 = DT.Omega(list(barrier1, 1:mesh1$t), mesh1)
# - creates a legal Omega object (subdivision of the mesh into areas)
Omega.SP1 = DT.polygon.omega(mesh1, Omega1)
# - creates polygons for the different areas: Barrier area and Normal area
# - the first polygon should be equal to the polygon p
# - it is important to check here that the polygons you created were actually the polygons you wanted
# -  - this cannot be taken for granted in the mesh ''extension''

## Visually check correctness
plot(mesh1)
plot(Omega.SP1[[1]], add=T, col='grey')
plot(Omega.SP1[[2]], add=T, col='lightblue')
plot(mesh1, add=T, draw.segments=F)
points(loc1)
# - check that the greay area represents the Barrier in a good way
# - and that the light blue area represents the normal terrain


### EXAMPLE 2: ELLIPSIS ISLAND ###
# You can think of this as an example of an island that is a barrier to spatial correlation of aquatic animals
# This example creates the mesh independently of the barrier area, and then adds the barrier afterwards
# - The barrier area is defined through distances, and not by polygons
# - This mesh is therefore NOT CONSISTENT when you decrease max.edge.length (i.e. it will look somewhat different)

## Mesh setup
loc2 = loc1
mesh2 = inla.mesh.2d(loc=loc2, max.e = max.edge.length, offset=1)

mesh2 = DT.mesh.addon.posTri(mesh2)
x = mesh2$posTri[ , 1]
# - is the x coordiate of the triangle positions
y = mesh2$posTri[ , 2]
# - is the y coordiate of the triangle positions

d = ((x-5)/2)^2 + ((y-5)/.5)^2
# - is the elliptical distance from an ellipse with centre in c(5,5) and scaling 2 and 0.5

barrier2 = which(d < 1^2)
# - is the triangle indices within the ellipse of scaled radius 1
Omega2 = DT.Omega(list(barrier2, 1:mesh2$t), mesh2)
# - creates a legal Omega object (subdivision of the mesh)
Omega.SP2 = DT.polygon.omega(mesh2, Omega2)
# - creates polygons for the different areas: Barrier area and Normal area

## Visually check correctness
plot(mesh2)
plot(Omega.SP2[[1]], add=T, col='grey')
plot(Omega.SP2[[2]], add=T, col='lightblue')
plot(mesh2, add=T)
points(loc2)
# - check that the black area represents the Barrier in a good way


### EXAMPLE 3: RIVER INLET ###
# This example is how I recommend that you create meshes
# It is slightly more involved
# This kind of mesh construction is also useful for the usual stationary models

## User input
max.edge.length = 0.04

## Extract a part of a map
# The following code is taken from example in ?map2SpatialPolygons
# You may use this package 'maps', or similar databases, or shapefiles, or other stuff
# As long as you are able to extract a SpatialPolygons variable describing the barrier

library('maps'); library('maptools')

if (TRUE) { # use high resolution map
  library('mapdata')
  data("worldHiresMapEnv")
  nor_coast_poly <- map("worldHires", "norway", fill=TRUE, col="transparent",
                        plot=FALSE, ylim=c(58,72))
} else { # use the low resolution map # if the other code fails to run
  nor_coast_poly <- map("world", "norway", fill=TRUE, col="transparent",
                        plot=FALSE, ylim=c(58,72))
  
}

# The following code is taken from example in ?map2SpatialPolygons
# - copy-paste begins
nor_coast_poly$names
IDs <- sapply(strsplit(nor_coast_poly$names, ":"), function(x) x[1])
nor_coast_poly_sp <- map2SpatialPolygons(nor_coast_poly, IDs=IDs,
                                         proj4string=CRS("+proj=longlat +datum=WGS84"))
# - copy-paste ended

## We are interested here in a small part of the Norwegian coast
# - it is important not to use a polygon for the whole Norwegian cost
# - limit the polygon to where you have data
# - this is to keep the mesh small
# - a too large mesh leads to a very slow algorithm!

loc3 = matrix(c(4,60.6, 8,60.6, 4,61.5, 8,61.5), 4, 2, byrow = T)
# - this is set to give the right extent of the mesh
# - you can use data locations instead
# - a bit extra offset will be added by the mesh function (later)
max.edge.length = 0.04
# - this is set to a length that is small compared to the width of the study area
square = bakka.square.polygon(xlim=range(loc3[ ,1]), ylim=range(loc3[ ,2]), 
                              ret.SP = T)
# - This is the study area
# - - you are interested in fitting and/or predicting locations inside here
# - You only care about the water area, since land is a barrier
# - If loc3 is your data locations, you may want to make this square slightly bigger
# - - You should never report plots that are bigger than this square
square@proj4string=nor_coast_poly_sp@proj4string
# - This creates the same coordinate system for the square as we have for the map
poly3 = gIntersection(nor_coast_poly_sp, square)
# - The part of Norway inside the study area square

plot(poly3, xlim=range(loc3[ ,1]), ylim=range(loc3[ ,2]), axes=T)
# - We see that poly3 is the land polygon
poly3.water = gDifference(square, poly3)
# - This is the polygon for water
# - This is what we need, since the mesh will be here!

mesh3 = inla.mesh.2d(boundary = poly3.water, max.e = max.edge.length)
# - A mesh over the water area, with the desired fidelity (max.e)
# - Many users fit (stationary) models to this mesh, but that is wrong!
# - - This is because you need some 'free space' between the study region and the 
#     border of the mesh for the finite element approximation to be correct
oldpar <- par(mar = c(1, 1, 3, 1))
plot(mesh3, main="Inner mesh")
mesh3$n
mesh3 = inla.mesh.2d(loc=mesh3$loc, max.edge = max.edge.length*2, offset=max.edge.length*2)
# - adds a convex hull with a small extension
# - here the mesh is coarser
# - this is done to have a good approximation on the border of the barrier
plot(mesh3, main="Convex mesh")
mesh3$n
mesh3 = inla.mesh.2d(loc=mesh3$loc, max.edge = max.edge.length*4, offset=0.5)
# - This adds an outer extension
# - This is needed since the spatial field is incorrect at the outer regions of
#   the mesh (inflated range and variance)
# - This outer extension will not be used for fitting or prediction, it is just
#   there to fix SPDE boundary approximation issues
plot(mesh3, main="Mesh with outer extension")
mesh3$n


mesh3 = DT.mesh.addon.posTri(mesh3)

## Attempt 1: Define which triangles that belong to land
# - This attempt is going to fail, but it is important that you know why, and how to check
points = SpatialPoints(mesh3$posTri)
points@proj4string = poly3@proj4string
barrier3 = over(poly3, points, returnList=T)
barrier3 = unlist(barrier3)
Omega3 = DT.Omega(list(barrier3, 1:mesh3$t), mesh3)
Omega.SP3 = DT.polygon.omega(mesh3, Omega3)

## Attempt 1: Visually check correctness
plot(mesh3)
plot(Omega.SP3[[1]], add=T, col='grey')
plot(Omega.SP3[[2]], add=T, col='lightblue')
plot(mesh3, add=T)
points(loc3, col='black')
# - Check that the grey area represents the Barrier in a good way
# - It does not! 

## Attempt 2: Define land triangles by using the big polygon
# - Previously we attempted to only use the polygon in our study area
# - But we extended the mesh!
barrier3 = over(nor_coast_poly_sp, points, returnList=T)
barrier3 = unlist(barrier3)
Omega3 = DT.Omega(list(barrier3, 1:mesh3$t), mesh3)
Omega.SP3 = DT.polygon.omega(mesh3, Omega3)

## Attempt 2: Visually check correctness
plot(mesh3)
plot(Omega.SP3[[1]], add=T, col='grey')
plot(Omega.SP3[[2]], add=T, col='lightblue')
plot(mesh3, add=T)
points(loc3, col='black')
# - check that the grey area represents the Barrier in a good way
# - It does!

## Comment: What if I do not have the larger polygon?
# - Then you can just draw some simple polygons to represent what you think the nearby coast is
# - For example, in this case, just continue the barrier vertically


### PROBLEM ###
# To check your understanding, figure out what is slightly wrong with the first mesh (the gap in a barrier), and change slightly the definition of the barrier so that it is more sensible.


### DEFINE local functions ###

## How to plot fields in 2D
local.plot.field = function(field, xlim, ylim , mesh, ...){
  proj = inla.mesh.projector(mesh, xlim = xlim, ylim = ylim, dims=c(300, 300))
  field.proj = inla.mesh.project(proj, field)
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), xlim = xlim, ylim = ylim, ...)  
}
local.plot.field.cutoff = function(field, xlim, ylim , mesh, ...){
  proj = inla.mesh.projector(mesh, xlim = xlim, ylim = ylim, dims=c(300, 300))
  field.proj = inla.mesh.project(proj, field)
  field.proj[field.proj<0.1] = NA
  # - E.g. correlation, we only plot when >= 0.1
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), xlim = xlim, ylim = ylim, ...)  
}

## How to compute the spatial correlation with respect to a location
# - Uses the precision matrix Q
local.find.correlation = function(Q, location, mesh) {
  d = sqrt(diag(inla.qinv(Q)))
  # - marginal standard deviations
  
  node.id = dt.nearest.mesh.index(location, mesh)
  node.loc = c(mesh$loc[node.id, ][1], mesh$loc[node.id, ][2])
  print(paste('The location used was c(', node.loc[1], ', ', node.loc[2], ')'))
  Inode = rep(0, dim(Q)[1]); Inode[node.id] = 1
  corr = solve(Q, Inode); 
  corr = drop(matrix(corr)) / d; corr = corr/d[node.id]
  
  return(corr)
}





### Heavy Computations ###
# This is where we solve the SPDEs
# This solution is dependent on the mesh construction
# The mesh construction answers 'How do you want to approximate the SPDE?'

print('Doing heavy pre-computations...')
spde.DT1 = DT.FEMmatrices(mesh1, Omega1)
spde.DT2 = DT.FEMmatrices(mesh2, Omega2)
spde.DT3 = DT.FEMmatrices(mesh3, Omega3)
print('Done')





### Example 1: Correlation - Illustrate in two locations ###
Q = DT.precision.new(spde.DT1, ranges = c(0.3, 3))

## First location
c = local.find.correlation(Q, location = c(5, 4.3), mesh=mesh1)
local.plot.field.cutoff(c, xlim = c(2, 8), ylim = c(2, 8), mesh=mesh1)
# - Plot the nicely colored correlation surface and the legend
plot(Omega.SP1[[1]], add=T, col='grey')
# - Plot the Barrier on top of this
# - How the field acts below this barrier is of no practical significance (you are not allowed to have data there)
# - - But you may find it academically interesting
points(mesh1$loc[which.max(c), c(1,2), drop=F], pch=19)

## Second location
c = local.find.correlation(Q, location = c(5, 4.9), mesh=mesh1)
local.plot.field.cutoff(c, xlim = c(2, 8), ylim = c(2, 8), mesh=mesh1)
plot(Omega.SP1[[1]], add=T, col='grey')
points(mesh1$loc[which.max(c), c(1,2), drop=F], pch=19)

## Plot the marginal standard deviation
# - It does not vary too much
d = sqrt(diag(inla.qinv(Q)))
local.plot.field(d, xlim=c(2, 8), ylim=c(2, 8), mesh=mesh1 )
plot(Omega.SP1[[1]], add=T, col='grey')

### Example 2: Correlation - Illustrate in two locations ###
Q = DT.precision.new(spde.DT2, ranges = c(0.3, 3))

## First location
c = local.find.correlation(Q, location = c(4, 4.7), mesh=mesh2)
local.plot.field.cutoff(c, xlim = c(2, 8), ylim = c(2, 8), mesh=mesh2)
plot(Omega.SP2[[1]], add=T, col='grey')
points(mesh2$loc[which.max(c), c(1,2), drop=F], pch=19)

## Second location
c = local.find.correlation(Q, location = c(7, 5), mesh=mesh2)
local.plot.field.cutoff(c, xlim = c(2, 8), ylim = c(2, 8), mesh=mesh2)
plot(Omega.SP2[[1]], add=T, col='grey')
points(mesh2$loc[which.max(c), c(1,2), drop=F], pch=19)

## Plot the marginal standard deviation
# - It does not vary too much
d = sqrt(diag(inla.qinv(Q)))
local.plot.field(d, xlim=c(2, 8), ylim=c(2, 8), mesh=mesh2 )
plot(Omega.SP2[[1]], add=T, col='grey')


### Example 3: Correlation - Illustrate in two locations ###
Q = DT.precision.new(spde.DT3, c(0.01, 0.5))

## First location
c = local.find.correlation(Q, loc=c(5.1, 61.1), mesh=mesh3)
local.plot.field.cutoff(c, xlim=range(loc3[ ,1]), ylim=range(loc3[ ,2]), mesh=mesh3)
plot(Omega.SP3[[1]], add=T, col='grey')
points(mesh3$loc[which.max(c), c(1,2), drop=F], pch=19)

## Second location
c = local.find.correlation(Q, loc=c(7, 61.1), mesh=mesh3)
local.plot.field.cutoff(c, xlim=range(loc3[ ,1]), ylim=range(loc3[ ,2]), mesh=mesh3)
plot(Omega.SP3[[1]], add=T, col='grey')
points(mesh3$loc[which.max(c), c(1,2), drop=F], pch=19)

## Plot the marginal standard deviation
# - This one varies quite a lot
# - But in a sensible way from the application point of view
# - If you do not want this, you can just rescale the Q-matrix, so that this sd is 1
d = sqrt(diag(inla.qinv(Q)))
local.plot.field(d, xlim=range(loc3[ ,1]), ylim=range(loc3[ ,2]), mesh=mesh3)
plot(Omega.SP3[[1]], add=T, col='grey')





### END ###
