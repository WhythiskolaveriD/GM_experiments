
### DESCRIPTION ###
# This code gives a minimal simulation-inference example for the Barrier model using R-INLA
# See: www.r-inla.org/barrier-model
# This example is only using simulated data, but you should be able to re-use parts of this code to run real datasets

### MASTER ###
# Ignore this! Internal use only!
# This document is either a master or a slave
# Master:
# - dtrcode repository
# Slaves: 
# - Gdrive/public/barriermodel1/Case162_SimulationInference_S.R
# Synchronise:
# - phdbakka/TutorialDT/TutorialDT.Rnw (parts only)


### LOADING LIBRARIES AND FUNCTIONS ###
rm(list=ls())
source('functMinimalDT.R')

### USER INPUT ###
max.edge.length = 0.4
# - The coarseness of the finite element approximation
# - Corresponds to grid-square width in discretisation/rasters
# - - Except that finite element approximations are much better
# - Should be compared to size of study area
# - Should be less than a third of the estimated spatial range
# - Computational time roughly *8 when you halve this value
# - - Except, in some cases, this *8 can be as low as *2
smalldist = 0.2
# - the width of the opening in the barrier
set.seed(2016)
set.inla.seed = 2016

### MESH: NARROW PASSAGE ###
# This code is the same as example 1 in Case161_CorrelationExample
# - EXAMPLE 1: NARROW PASSAGE

width = .5
# - The width/thickness of the barrier
p2 = bakka.square.polygon(xlim=c(-1, 5-smalldist/2), ylim=5+width*c(-.5, .5), 
                          ret.SP = T)
p3 = bakka.square.polygon(xlim=c(5+smalldist/2, 11), ylim=5+width*c(-.5, .5), 
                          ret.SP = T)
p = SpatialPolygons(c(p2@polygons, p3@polygons))
# - The total barrier area polygon
# - Use plot(p) to investigate

loc1 = matrix(c(0,0, 10,0, 0,10, 10,10), 4, 2, byrow = T)
# - This defines the extent of the mesh
# - In an application, if you want the mesh to depend on your data locations, 
#   you may use those locations here
seg = inla.sp2segment(p)
mesh = inla.mesh.2d(loc=loc1, interior = seg, max.e = max.edge.length, 
                    offset=1)
# - The INLA mesh constructor, used for any INLA-SPDE model
mesh = DT.mesh.addon.posTri(mesh)
# - A special tool needed for these models

barrier = over(p, SpatialPoints(mesh$posTri), returnList=T)
# - checking which mesh triangles are inside the barrier area
barrier = unlist(barrier)
Omega = DT.Omega(list(barrier, 1:mesh$t), mesh)
Omega.SP = DT.polygon.omega(mesh, Omega)
# - creates polygons for the different areas: Barrier area and Normal area

## Visually check correctness
oldpar <- par(mar = c(1,1,3,1))
plot(mesh, main="")
plot(Omega.SP[[1]], add=T, col='grey')
plot(Omega.SP[[2]], add=T, col='lightblue')
plot(mesh, add=T)
points(loc1)






### DEFINE how to plot the field ###
local.plot.field = function(field, ...){
  xlim = c(2, 8); ylim = xlim;
  proj = inla.mesh.projector(mesh, xlim = xlim, ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300 by 300 grid for plotting
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, ...)  
  # - Use image.plot to get nice colors and legend
}
print(mesh$n)
# - This is the appropriate length of the field variable

### Heavy Computations ###
print('Doing heavy pre-computations...')
spde.DT = DT.FEMmatrices(mesh, Omega)
# - Set up the matrices for solving the SPDE
# - May be time consuming
print('Done')









### SIMULATE DATA ###

## User input
ranges = c(0.3, 3)
# - the first range is for the barrier area
# - - it is not sensitive to the exact value here, just make it ''small''
# - the second range is for the normal area
sigma.u = 1
# - size of the random effect
sigma.epsilon = 0.2
# - size of the iid noise in the Gaussian likelihood

## Sample the true field
Q = DT.precision.new(spde.DT, ranges)
# - Construct the precision matrix for simulating

u = inla.qsample(n=1, Q=Q, seed = set.inla.seed)
# - simulate a spatial field u with this precision matrix
u = u[ ,1]

local.plot.field(u)
plot(Omega.SP[[1]], add=T, col='grey')
title(main='True spatial field')

## Sample the data locations
num.try = 500 
# - try to sample this number of data locations
# - locations sampled inside the barrier will be removed in a few lines
loc.data.before = matrix(runif(num.try*2, min=2, max=8), num.try, 2)
temp = SpatialPoints(loc.data.before)
ok.locations = !is.na(over(temp, Omega.SP[[2]]))
# - only allow locations that are not inside the Barrier area
loc.data = loc.data.before[ok.locations, ]
names(loc.data) = c('locx', 'locy')

A.data = inla.spde.make.A(mesh, loc.data)
# - the projector matrix required for any spatial model
u.data = A.data %*% u
# - projecting the field from the finite element representation to the data locations
df = data.frame(loc.data)
# - df is the dataframe used for modeling
df$y = drop(sigma.u*u.data + sigma.epsilon*rnorm(nrow(df)))
# - sample observations with gaussian noise

stk <- inla.stack(data=list(y=df$y), A=list(A.data, 1),
                  effects=list(s=1:mesh$n, intercept=rep(1, nrow(df))), 
                  remove.unused = FALSE, tag='est')
# - this is the common stack used in INLA SPDE models
# - see the SPDE-tutorial
# - - http://www.r-inla.org/examples/tutorials/spde-tutorial











### Stationary model ###
# - This is not the new model, but the common spatial model in INLA
# - We add this for comparison between stationary and Barrier model

spde.model = inla.spde2.matern(mesh)
# - Set up the model component for the spatial SPDE model: Stationary Matern model
# - I assume you are familiar with this model

formula <- y ~ 0 + intercept + f(s, model=spde.model)
# - Remove the default intercept
# - - Having it in the stack instead improves the numerical accuracy of the
#     INLA algorithm
# - Fixed effects + random effects

print('Running inference with a stationary model...')
res.stationary <- inla(formula, data=inla.stack.data(stk),
            control.predictor=list(A = inla.stack.A(stk)),
            family = 'gaussian',
            control.family = list(hyper = list(prec = list(
              prior = "pc.prec", fixed = FALSE, param = c(0.2,0.5)))),
            control.compute = list(config=T, dic=T, cpo=T, waic=T),
            control.mode=list(restart=T, theta=c(4,-1.7,0.25))
            # - Initial values for the hyperparameters for the numerical optimisation
            # - Internal scale: log(precision, log(sd), log(range)
            # - Starting with good values makes the algorithm faster
            # - You can remove this
            # - If you pick very bad values, the numerics in inla(...) can crash
            )
print('Done')

## Investigate results
summary(res.stationary)

local.plot.field(res.stationary$summary.random$s$mean)
# - plot the posterior spatial marginal means
# - we call this the spatial estimate, or the smoothed data
plot(Omega.SP[[1]], add=T, col='grey')
# - Posterior spatial estimate using the stationary model
title(main='Posterior spatial estimate using the stationary model')








#TODO from here
### Barrier model ###
# This is the new model!
# The Barrier model is a special case of the Different Terrains model, therefore the naming is DT in variables and functions

DT = DT.rgeneric.environment(spde.DT, prior.parameters=c(1, 1), 
                             fixed.ranges = c(0.05, NA))
# - This function creates the environment/variables needed to run the internal 
#   rgeneric engine inside the C-code in R-INLA
# - We fix the barrier range to a different value than we used for simulations
# - - Why? It does not matter, as long as it is 'small' the models are very
#     similar
# - - This shows that you do not need to know the true 'barrier range'!
# - The prior parameters are the lambdas in the exponential priors for standard 
#   deviation and inverse-range

filename.rgeneric.environment = 'temp.rgeneric.init.RData'
save('DT', file = filename.rgeneric.environment)
# - save the new environment file
# - this is not an actual R environment, just a list-of-lists that is loaded 
#   at rgeneric startup

barrier.model = inla.rgeneric.define(DT.rgeneric.dt.model, 
                                     n=dim(spde.DT$I)[1], ntheta = DT$ntheta, 
                                     R.init=filename.rgeneric.environment, debug=F) 
# - This contains all the needed functions for the spatial model component
# - - Most important is the function computing the precision matrix Q from the 
#     hyperparameters theta
# - These functions are called by the C-code in R-INLA when R-INLA is 
#   performing inference
# - - The inla(...) call spawns a C algorithm running inference, which again
#     spawns an R algorithm running rgeneric


formula2 <- y ~ 0 + intercept + f(s, model=barrier.model)
# - The spatial model component is different from before
# - The rest of the model setup is the same! (as in the stationary case)
# - - e.g. the inla(...) call below is the same, only this formula is different

print('Running inference with the new Barrier model...')
res.barrier <- inla(formula2, data=inla.stack.data(stk),
                       control.predictor=list(A = inla.stack.A(stk)),
                       family = 'gaussian',
                       control.family = list(hyper = list(prec = list(
                         prior = "pc.prec", fixed = FALSE, param = c(0.2,0.5)))),
                       control.compute = list(config=T, dic=T, cpo=T, waic=T),
                    control.mode=list(restart=T, theta=c(3.2, 0.4, 1.6)),
                    verbose=T)
print('Done')

## Investigate results 
summary(res.barrier)

local.plot.field(res.barrier$summary.random$s$mean)
# - plot the posterior spatial marginal means
# - we call this the spatial (smoothing) estimate
plot(Omega.SP[[1]], add=T, col='grey')
title(main='Posterior spatial estimate using the Barrier model')




### MODEL COMPARISON ###
dic1 = c(res.stationary$dic$dic, res.barrier$dic$dic)
waic1 = c(res.stationary$waic$waic, res.barrier$waic$waic)
model.comparison = rbind(dic1, waic1)
rownames(model.comparison) = c('dic', 'waic')
colnames(model.comparison) = c('stationary', 'barrier')
print(model.comparison)
# - We see that the barrier model is better in both cases
# - - With the original seeds
# - - Different seeds can give different results
# - There is a large amount of ''randomness'' in these numbers,
#   so please do not take them too seriously



### END ###
