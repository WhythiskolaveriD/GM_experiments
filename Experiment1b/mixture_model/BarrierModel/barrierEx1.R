
## ---- warning=FALSE, message=FALSE---------------------------------------
library(INLA); library(sp); library(fields)

## Load data
#load(file = "https://haakonbakka.bitbucket.io/data/WebSiteData-Archipelago.RData")
# - if you have an internet connection
load(file = "WebSiteData-Archipelago.RData")
# - if you have saved the file locally

## What is loaded
# - poly.water is our study area
# - df is our dataframe to be analysed
# - dat is the orginial dataframe

## Load functions
#source('https://haakonbakka.bitbucket.io/functions-barriers-dt-models-march2017.R')
# - if you have an internet connection
source('functions-barriers-dt-models-march2017.R')
# - if you have saved the file locally

set.seed(2016)
set.inla.seed = 2016

## ------------------------------------------------------------------------
max.edge = 0.6
bound.outer = 4.6
mesh = inla.mesh.2d(boundary = poly.water,
                    loc=cbind(df$locx, df$locy),
                    max.edge = c(1,5)*max.edge,
                    cutoff = 0.06,
                    offset = c(max.edge, bound.outer))
plot(mesh, main="Our mesh", lwd=0.5); points(df$locx, df$locy, col="red")
mesh$n

## ------------------------------------------------------------------------
df$y = df$y.smelt

## ------------------------------------------------------------------------
M = list()
# - the list of all models
M[[1]] = list()
# - the first model
M[[1]]$shortname = "stationary-no-cov"
# - a short description, the name of the model
M[[2]] = list()
M[[2]]$shortname = "stationary-all-cov"
M[[3]] = list()
M[[3]]$shortname = "barrier-no-cov"
M[[4]] = list()
M[[4]]$shortname = "barrier-all-cov"

## ------------------------------------------------------------------------
A.i.s = inla.spde.make.A(mesh, loc=cbind(df$locx, df$locy))

## ---- warning=FALSE, message=FALSE---------------------------------------
stk <- inla.stack(data=list(y=df$y, e=df$exposure), 
                  effects=list(s=1:mesh$n,
                               data.frame(m=1, df[ ,5:11]), 
                               # - m is the intercept
                               iidx=1:nrow(df)),
                  A=list(A.i.s, 1, 1),
                  remove.unused = FALSE, tag='est')  

## ---- warning=FALSE, message=FALSE---------------------------------------
spde = inla.spde2.pcmatern(mesh, prior.range = c(6, .5), prior.sigma = c(3, 0.01))
# - We put the prior median at approximately 0.5*diff(range(df$locy))
# - - this is roughly the extent of our study area
# - The prior probability of marginal standard deviation 3 or more is 0.01.

## ------------------------------------------------------------------------
hyper.iid = list(prec = list(prior = 'pc.prec', param = c(3, 0.01))) 
# - the param the same as prior.sigma above, with the same interpretation
M[[1]]$formula = y~ -1+m + f(s, model=spde) + f(iidx, model="iid", hyper=hyper.iid)
# - no covariates (except intercept m)
M[[2]]$formula = as.formula(paste( "y ~ -1 + ",paste(colnames(df)[5:11], collapse = " + ")))
M[[2]]$formula = update(M[[2]]$formula, .~. +m + f(s, model=spde) + f(iidx, model="iid", hyper=hyper.iid))

## ------------------------------------------------------------------------
print(M[[1]])
print(M[[2]])

## ------------------------------------------------------------------------
mesh = dt.mesh.addon.posTri(mesh)
# - compute the triangle positions
normal = over(poly.water, SpatialPoints(mesh$posTri), returnList=T)
# - checking which mesh triangles are inside the normal area
normal = unlist(normal)
Omega = dt.Omega(list(normal, 1:mesh$t), mesh)
Omega.SP = dt.polygon.omega(mesh, Omega)

## ------------------------------------------------------------------------
Q.barrier = dt.create.Q(mesh, Omega, copy.ranges = c(NA, 1),
                        copy.ranges.frac = c(1, 1/5))
# - We let the second range be a copy of the first range, multiplied by 1/5
# - Why? The value 1/5 is our default recommendation for the fidelity of the spatial approximation
# - Choosing a fraction other than 1/5 does not impact the results significantly. This shows that you do not need to know the true 'barrier range'
# - time: Ca 1 min

log.prior = dt.create.prior.log.exp(
  prior.param = c(-log(0.01)/3, -log(0.5)*6))
# - The prior parameters are the lambdas in the exponential 
#   priors for standard deviation and inverse-range
# - the first is log(prob)/exceed, the second log(prob)*exceed
# - the second is exponential for inverse range, therefore multiplication!

barrier.model = dt.inla.model(
  Q = Q.barrier, log.prior=log.prior)


## ------------------------------------------------------------------------
M[[3]]$formula = y~ -1+m + f(s, model=barrier.model) + f(iidx, model="iid", hyper=hyper.iid)
# - no covariates (except intercept m)
M[[4]]$formula = as.formula(paste( "y ~ -1 + ",paste(colnames(df)[5:11], collapse = " + ")))
M[[4]]$formula = update(M[[4]]$formula, .~. +m + f(s, model=barrier.model) + f(iidx, model="iid", hyper=hyper.iid))

## ------------------------------------------------------------------------
print(M[[3]])
print(M[[4]])

## ------------------------------------------------------------------------
## Initial values
# - speeds up computations
# - improves accuracy of computations
# - set these to NULL the first time you run the model
M[[1]]$init = c(2.509,1.199,-0.574)
M[[2]]$init = c(1.162,0.313,-0.627)
M[[3]]$init = c(0.833,2.244,-0.471)
M[[4]]$init = c(0.044,1.274,-0.596)

## ---- warning=FALSE, message=FALSE---------------------------------------
for (i in 1:length(M)){
  print(paste("Running:  ", M[[i]]$shortname))
  M[[i]]$res = inla(M[[i]]$formula,
                    data=inla.stack.data(stk),
                    control.predictor=list(A=inla.stack.A(stk)),
                    family="poisson", E = e,
                    control.inla= list(int.strategy = "eb"),
                    control.mode=list(restart=T, theta=M[[i]]$init), verbose = TRUE)  
}
# - time: < 20 min

## ------------------------------------------------------------------------
for (i in 1:length(M)){
  print(paste(round(M[[i]]$res$internal.summary.hyperpar$mode, 3), collapse = ','))
}

## ------------------------------------------------------------------------
summary(M[[1]]$res)

## ------------------------------------------------------------------------
summary(M[[2]]$res)

## ------------------------------------------------------------------------
summary(M[[3]]$res)

## ------------------------------------------------------------------------
summary(M[[4]]$res)

## ------------------------------------------------------------------------
#M[[i]]$res$logfile

## ------------------------------------------------------------------------
local.plot.field = function(field, mesh, xlim, ylim, ...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  if (missing(xlim)) xlim = poly.water@bbox[1, ] 
  if (missing(ylim)) ylim = poly.water@bbox[2, ]
  # - choose plotting region to be the same as the study area polygon
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, ...)  
}

## ------------------------------------------------------------------------
for (i in c(1,3)) {
  field = M[[i]]$res$summary.random$s$mean + M[[i]]$res$summary.fixed['m', 'mean']
  local.plot.field(field, mesh, main=paste(M[[i]]$shortname), zlim=c(-10.5, 1))
  plot(Omega.SP[[2]], add=T, col="grey")
  points(df$locx, df$locy)
}

## ------------------------------------------------------------------------
for (i in c(2,4)) {
  field = M[[i]]$res$summary.random$s$mean + M[[i]]$res$summary.fixed['m', 'mean']
  local.plot.field(field, mesh, main=paste(M[[i]]$shortname), zlim=c(-9, -4.5))
  plot(Omega.SP[[2]], add=T, col="grey")
  points(df$locx, df$locy)
}

## ------------------------------------------------------------------------
## Set up dataframe with relevant results
res = M[[2]]$res$summary.fixed[ ,c(4,3,5)]
# - all covar, some quantiles
for(i in c(4)) res = rbind(res, M[[i]]$res$summary.fixed[ ,c(4,3,5)])

colnames(res) = c("E", "L", "U")
rownames(res)=NULL
n.covar = nrow(M[[2]]$res$summary.fixed)
res$model = factor(rep(c("MS", "MB"), each=n.covar, 
                       levels = c("MS", "MB")))
res$covar = factor(rep(rownames(M[[2]]$res$summary.fixed), 2))

if(require(ggplot2)) {
  ggplot(res, aes(x = model, y = E)) +
    facet_wrap(~covar, scales = "free_y") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymax = U, ymin = L)) +
    xlab(NULL) + ylab(NULL) 
}
