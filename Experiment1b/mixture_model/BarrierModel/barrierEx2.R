

## ---- warning=FALSE, message=FALSE---------------------------------------
library(INLA)
library(mgcv)

#source('https://haakonbakka.bitbucket.io/functions-barriers-dt-models-march2017.R')
# - if you have an internet connection
source('functions-barriers-dt-models-march2017.R')
# - if you have saved the file locally

set.seed(2016)
set.inla.seed = 2016

## ------------------------------------------------------------------------
N.loc = 100
# - number of locations
# - 100 in the soap-film paper (Wood)
sigma.eps = 0.1
# - measurement noise
# - Wood uses 0.1, 1, and 10
global.zlim=c(-1, 1)*4.3

## ------------------------------------------------------------------------
## plot the function, and its boundary...
fsb <- fs.boundary()
m<-300;n<-150 
xm <- seq(-1,4,length=m);yn<-seq(-1,1,length=n)
xx <- rep(xm,n);yy<-rep(yn,rep(m,n))
tru = fs.test(xx,yy, b = 1)
tru.matrix <- matrix(tru,m,n) ## truth
image.plot(xm,yn,tru.matrix,xlab="x",ylab="y", asp=1)
lines(fsb$x,fsb$y,lwd=3)
contour(xm,yn,tru.matrix,levels=seq(global.zlim[1], global.zlim[2],len=20),add=TRUE, col="white", drawlabels = F)
range(tru, na.rm = T)

## ------------------------------------------------------------------------
dat = data.frame(y = tru, locx = xx, locy=yy)
dat = dat[!is.na(dat$y), ]
df = dat[sample(1:nrow(dat), size=100) ,]
df$y = df$y + rnorm(100)*sigma.eps
str(df)
summary(df)

## ------------------------------------------------------------------------
p = Polygon(cbind(fsb$x, fsb$y))
p = Polygons(list(p), ID = "none")
poly = SpatialPolygons(list(p))
plot(poly)
points(df$locx, df$locy)
axis(1); axis(2)

## ------------------------------------------------------------------------
max.edge = 0.2
bound.outer = 1.5
mesh = inla.mesh.2d(boundary = poly,
                    loc=cbind(df$locx, df$locy),
                    max.edge = c(1,5)*max.edge,
                    #cutoff = 0.1,
                    cutoff = 0.04,
                    # 0.1 is fast and bad, 0.04 ok?
                    offset = c(max.edge, bound.outer))

plot(mesh, main="Our mesh", lwd=0.5)
mesh$n

## ------------------------------------------------------------------------
A.i.s = inla.spde.make.A(mesh, loc=cbind(df$locx, df$locy))
stk = inla.stack(data=list(y=df$y), 
                 effects=list(s=1:mesh$n,
                              m = rep(1, nrow(df))),
                 A=list(A.i.s, 1),
                 remove.unused = FALSE, tag='est') 

## ---- warning=FALSE, message=FALSE---------------------------------------
prior.range = c(1, .5)
prior.sigma = c(3, 0.01)
spde = inla.spde2.pcmatern(mesh, prior.range=prior.range, prior.sigma=prior.sigma)
# - We put the prior median at approximately 0.5*diff(range(df$locy))
# - - this is roughly the extent of our study area
# - The prior probability of marginal standard deviation 3 or more is 0.01.

## ------------------------------------------------------------------------
M = list()
M[[1]] = list(shortname="stationary-model")
M[[1]]$formula = y~ -1+m + f(s, model=spde)

## ------------------------------------------------------------------------
mesh = dt.mesh.addon.posTri(mesh)
# - compute the triangle positions
normal = over(poly, SpatialPoints(mesh$posTri), returnList=T)
# - checking which mesh triangles are inside the normal area
normal = unlist(normal)
Omega = dt.Omega(list(normal, 1:mesh$t), mesh)
Omega.SP = dt.polygon.omega(mesh, Omega)
plot(Omega.SP[[2]], col="grey", main="The barrier region (in grey)")

## ------------------------------------------------------------------------
Q.barrier = dt.create.Q(mesh, Omega, copy.ranges = c(NA, 1),
                        copy.ranges.frac = c(1, 1/5))
# - [TODO]
# - 1/5 of the range is our general recommendation for "resolution" of the "spatial computations"
# - - e.g. maximum recommended max.edge.length
# - time: Ca 10 sec

log.prior = dt.create.prior.log.exp(
  prior.param = c(-log(prior.sigma[2])/prior.sigma[1], -log(prior.range[2])/prior.range[1])) 
#c(-log(0.01)/3, -log(0.5)*6))
# - The prior parameters are the lambdas in the exponential 
#   priors for standard deviation and inverse-range
# - the first is log(prob)/exceed, the second log(prob)*exceed
# - the second is exponential for inverse range, therefore multiplication!

barrier.model = dt.inla.model(
  Q = Q.barrier, log.prior=log.prior)


## ------------------------------------------------------------------------
M[[2]] = list(shortname="barrier-model")
M[[2]]$formula = y~ -1+m + f(s, model=barrier.model)

## ------------------------------------------------------------------------
## Initial values
# - speeds up computations
# - improves accuracy of computations
# - set these to NULL the first time you run the model
M[[1]]$init = c(4.942,1.006,0.851)
M[[2]]$init = c(4.511,0.521,2.391)

## ---- warning=FALSE, message=FALSE---------------------------------------
hyper.iid = list(prec = list(prior = 'pc.prec', param = c(3, 0.01))) 
hyper.fixed = list(prec = list(init=2.3, fixed=T))
# - sig = exp(-0.5*logprec)
# - sig=0.01, logsig=-4.6, logprec = 2.3

start.time <- Sys.time()
for (i in 1:length(M)){
  print(paste("Running:  ", M[[i]]$shortname))
  M[[i]]$res = inla(M[[i]]$formula,
                    data=inla.stack.data(stk),
                    control.predictor=list(A=inla.stack.A(stk)),
                    family="gaussian", 
                    control.family = list(hyper = hyper.iid),
                    #control.family = list(hyper = hyper.fixed),
                    control.inla= list(int.strategy = "eb"),
                    #verbose=T,
                    control.mode=list(restart=T, theta=M[[i]]$init))  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
# - time: ca 30 sec

## ------------------------------------------------------------------------
for (i in 1:length(M)){
  print(paste(round(M[[i]]$res$internal.summary.hyperpar$mode, 3), collapse = ','))
}

## ------------------------------------------------------------------------
summary(M[[1]]$res)

## ------------------------------------------------------------------------
summary(M[[2]]$res)

## ------------------------------------------------------------------------
#M[[i]]$res$logfile

## ------------------------------------------------------------------------
local.plot.field = function(field, mesh, xlim, ylim, zlim, n.contours=10, ...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  if (missing(xlim)) xlim = poly@bbox[1, ] 
  if (missing(ylim)) ylim = poly@bbox[2, ]
  # - choose plotting region to be the same as the study area polygon
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  if (missing(zlim)) zlim = range(field.proj)
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, asp=1, ...)  
  contour(x = proj$x, y=proj$y, z = field.proj,levels=seq(zlim[1], zlim[2],length.out = n.contours),add=TRUE, drawlabels=F, col="white")
  # - without contours it is very very hard to see what are equidistant values
}

## ------------------------------------------------------------------------
for (i in c(1,2)) {
  field = M[[i]]$res$summary.random$s$mean + M[[i]]$res$summary.fixed['m', 'mean']
  local.plot.field(field, mesh, main=paste(M[[i]]$shortname), zlim=global.zlim, n.contours = 20)
  plot(Omega.SP[[2]], add=T, border="black", col="white")
  points(df$locx, df$locy)
}

## ------------------------------------------------------------------------
## Truth on the grid
summary(dat)

## Remember
# M[[1]] is the stationary, M[[2]] is the barrier model

A.grid = inla.spde.make.A(mesh, loc=cbind(dat$locx, dat$locy))
for (i in 1:2) {
  M[[i]]$est = drop(A.grid %*% M[[i]]$res$summary.random$s$mean) +
    M[[i]]$res$summary.fixed["m", "mean"]
  M[[i]]$rmse = sqrt(mean((M[[i]]$est-dat$y)^2))
  M[[i]]$mae = mean(abs(M[[i]]$est-dat$y))
  M[[i]]$mae.sd = sd(abs(M[[i]]$est-dat$y))/sqrt(length(M[[i]]$est))
}

## Display results
data.frame(name=unlist(lapply(M, function(x) c(x$shortname))),
           rmse=unlist(lapply(M, function(x) c(x$rmse))),
           mae=unlist(lapply(M, function(x) c(x$mae))),
           mae.sd=unlist(lapply(M, function(x) c(x$mae.sd))))
