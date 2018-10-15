### R code from vignette source 'interface.Rtex'

###################################################
### code chunk number 1: interface.Rtex:1-6
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("INLA")
library("lattice")
set.seed(12345L)
inla.qsample(n = 1, Matrix(1, 1, 1), seed = 12345L)


###################################################
### code chunk number 2: interface.Rtex:17-43
###################################################
my.levelplot <-
  function(proj, values,
           col.regions = grey.colors(64, 1, 0),
           aspect = "fill",
           contour = TRUE, labels = FALSE,
           xlim = range(proj$x), ylim = range(proj$y),
           ...) {
    z = inla.mesh.project(proj, values)
    print(levelplot(
      row.values = proj$x, column.values = proj$y, x = z,
      xlim = xlim, ylim = ylim,
      col.regions = col.regions, aspect = aspect,
      contour = contour, labels = labels, ...))
}
lattice.to.plot <- function(latt, projection, skip = 10) {
  loc = inla.mesh.map(latt$loc, projection = projection, inverse = FALSE)
  mat1 = matrix(loc[, 1], nrow = latt$dims[1], ncol = latt$dims[2])
  mat2 = matrix(loc[, 2], nrow = latt$dims[1], ncol = latt$dims[2])
  hskip = seq(1, latt$dims[1], skip)
  vskip = seq(1, latt$dims[2], skip)
  return(rbind(
    cbind(as.vector(rbind(cbind(mat1[, vskip], NA), NA)),
          as.vector(rbind(cbind(mat2[, vskip], NA), NA))),
    cbind(as.vector(rbind(cbind(t(mat1[hskip, ]), NA), NA)),
          as.vector(rbind(cbind(t(mat2[hskip, ]), NA), NA)))))
}


###################################################
### code chunk number 3: interface.Rtex:57-59 (eval = FALSE)
###################################################
## source("http://www.math.ntnu.no/inla/givemeINLA-testing.R")
## library("INLA")


###################################################
### code chunk number 4: interface.Rtex:63-64 (eval = FALSE)
###################################################
## inla.upgrade(testing = TRUE)


###################################################
### code chunk number 5: interface.Rtex:93-108
###################################################
mesh1d.basis.plot = function(mesh, n, idx = 1:mesh$m, w = NULL, add = FALSE, ...) {
  x = seq(mesh$interval[1], mesh$interval[2], length = n)
  A = inla.mesh.1d.A(mesh, x)
  if (is.null(w)) {
    if (!add)
      plot(x = 0, type = "n", xlim = mesh$interval, ylim = range(A), ...)
    for (k in idx)
      lines(x, A[,k])
  } else {
    f = A %*% w
    if (!add)
      plot(x = 0,type = "n", xlim = mesh$interval, ylim = range(f), ...)
    lines(x, f)
  }
}


###################################################
### code chunk number 6: interface.Rtex:121-132 (eval = FALSE)
###################################################
## inla.mesh.create(loc=, ...)
## inla.mesh.create(loc=, tv=, ...)
## inla.mesh.create(lattice=, ...)
## inla.mesh.create(globe=, ...)
## inla.delaunay(loc=)
## inla.mesh.2d(loc=, ...)
## inla.mesh.1d(loc=, ...)
## 
## inla.sp2segment(sp)
## inla.mesh.segment(loc, ...)
## inla.nonconvex.hull(loc, ...)


###################################################
### code chunk number 7: interface.Rtex:156-160
###################################################
m = 50
points = matrix(runif(m * 2), m, 2)
mesh = inla.mesh.2d(
  loc = points, cutoff = 0.05, offset = c(0.1, 0.4), max.edge = c(0.05, 0.5) )


###################################################
### code chunk number 8: mesh1
###################################################
plot(mesh, main = "")
points(points)


###################################################
### code chunk number 9: interface.Rtex:186-187
###################################################
bnd = inla.nonconvex.hull(points, convex = 0.12)


###################################################
### code chunk number 10: mesh2
###################################################
lines(bnd, add = FALSE)
points(points)


###################################################
### code chunk number 11: interface.Rtex:195-196
###################################################
mesh = inla.mesh.2d(boundary = bnd, cutoff = 0.05, max.edge = c(0.1) )


###################################################
### code chunk number 12: mesh3
###################################################
plot(mesh, main = "")


###################################################
### code chunk number 13: interface.Rtex:202-204
###################################################
mesh = inla.mesh.2d(
  boundary = bnd, cutoff = 0.05, offset = c(1, 0.5), max.edge = c(0.1, 0.5) )


###################################################
### code chunk number 14: mesh4
###################################################
plot(mesh, main = "")


###################################################
### code chunk number 15: interface.Rtex:239-241 (eval = FALSE)
###################################################
## loc.cartesian = inla.mesh.map(loc.longlat, projection = "longlat")
## mesh2 = inla.mesh.2d(loc = loc.cartesian, ...)


###################################################
### code chunk number 16: interface.Rtex:244-245
###################################################
mesh2 = inla.mesh.create(globe = 10)


###################################################
### code chunk number 17: interface.Rtex:260-261 (eval = FALSE)
###################################################
## A = inla.spde.make.A(mesh, loc = points)


###################################################
### code chunk number 18: interface.Rtex:286-287
###################################################
spde = inla.spde2.matern(mesh, alpha = 2)


###################################################
### code chunk number 19: interface.Rtex:300-310
###################################################
sigma0 = 1
size = min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
range0 = size / 5
kappa0 = sqrt(8) / range0
tau0 = 1 / (sqrt(4 * pi) * kappa0 * sigma0)
spde = inla.spde2.matern(mesh,
  B.tau = cbind(log(tau0), -1, +1),
  B.kappa = cbind(log(kappa0), 0, -1),
  theta.prior.mean = c(0, 0),
  theta.prior.prec = c(0.1, 1) )


###################################################
### code chunk number 20: interface.Rtex:358-360
###################################################
truevar = (3 * sigma0)^2
truerange = range0


###################################################
### code chunk number 21: interface.Rtex:362-363
###################################################
Q = inla.spde.precision(spde, theta=c(log(3), 0))


###################################################
### code chunk number 22: interface.Rtex:366-367 (eval = FALSE)
###################################################
## x = inla.qsample(n = 2, Q)


###################################################
### code chunk number 23: interface.Rtex:370-371
###################################################
x = inla.qsample(n = 2, Q, seed = 123L)


###################################################
### code chunk number 24: interface.Rtex:383-384 (eval = FALSE)
###################################################
## x = inla.qsample(n = 2, Q, constr = spde$f$extraconstr)


###################################################
### code chunk number 25: interface.Rtex:401-404
###################################################
plot(mesh)
plot(mesh, rgl = TRUE)
lines(mesh$segm$bnd, mesh$loc, add = FALSE)


###################################################
### code chunk number 26: interface.Rtex:413-416
###################################################
plot(mesh, rgl = TRUE, col = x[, 1],
     color.palette = function(n) grey.colors(n, 1, 0),
     draw.edges = FALSE, draw.segments = TRUE, draw.vertices = FALSE)


###################################################
### code chunk number 27: interface.Rtex:422-424
###################################################
proj = inla.mesh.projector(mesh, dims = c(100, 100))
image(proj$x, proj$y, inla.mesh.project(proj, field = x[, 1]))


###################################################
### code chunk number 28: interface.Rtex:435-438
###################################################
mesh2 = inla.mesh.create(globe = 10)
proj2a = inla.mesh.projector(mesh2, projection = "longlat", dims = c(361, 181))
proj2b = inla.mesh.projector(mesh2, projection = "mollweide", dims = c(361, 181))


###################################################
### code chunk number 29: interface.Rtex:442-445
###################################################
spde2 = inla.spde2.matern(mesh2)
Q2 = inla.spde2.precision(spde2, theta = c(0, -1))
x2 = inla.qsample(n = 1, Q2, seed = 1234L)[, 1]


###################################################
### code chunk number 30: sample1
###################################################
my.levelplot(proj, x[, 1], at = pretty(x, 16), aspect = "iso",
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "", ylab = "", main = "")


###################################################
### code chunk number 31: sample2
###################################################
my.levelplot(proj, x[, 2], at = pretty(x, 16), aspect = "iso",
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "", ylab = "", main = "")


###################################################
### code chunk number 32: interface.Rtex:471-475
###################################################
latt=inla.mesh.lattice(
  x = seq(-180, 180, 1),
  y = seq(-90, 90, 1) * (1 - 1e-8),
  units = "longlat")


###################################################
### code chunk number 33: interface.Rtex:479-487
###################################################
print(contourplot(
  x = matrix(c(0, 0, 0, 0), 2, 2),
  row.values = c(0, 1), column.values = c(0, 1),
  xlim = c(-180, 180), ylim = c(-90, 90),
  aspect = 1/2, xlab = "Longitude", ylab = "Latitude"))
trellis.focus("panel", 1, 1, highlight = FALSE)
print(llines(lattice.to.plot(latt, "longlat"), col = 1))
trellis.unfocus()


###################################################
### code chunk number 34: interface.Rtex:490-498
###################################################
print(contourplot(
  x = matrix(c(0, 0, 0, 0), 2, 2),
  row.values = c(0, 1), column.values = c(0, 1),
  xlim = c(-180, 180), ylim = c(-1, 1),
  aspect = 1/2, xlab = "Longitude", ylab = ""))
trellis.focus("panel", 1, 1, highlight = FALSE)
print(llines(lattice.to.plot(latt, "longsinlat"), col = 1))
trellis.unfocus()


###################################################
### code chunk number 35: interface.Rtex:501-509
###################################################
print(contourplot(
  x = matrix(c(0, 0, 0, 0), 2, 2),
  row.values = c(0, 1), column.values = c(0, 1),
  xlim = c(-2, 2), ylim=c(-1, 1),
  aspect = 1/2, xlab = "", ylab = ""))
trellis.focus("panel", 1, 1, highlight = FALSE)
print(llines(lattice.to.plot(latt, "mollweide"), col = 1))
trellis.unfocus()


###################################################
### code chunk number 36: proj2a
###################################################
my.levelplot(proj2a, x2, at = pretty(x2, 16), aspect = 1/2,
      xlab = "Longitude", ylab = "Latitude", main = "Lon-Lat projection")
trellis.focus("panel", 1, 1, highlight = FALSE)
print(llines(lattice.to.plot(latt, "longlat", 30), col = 1))
trellis.unfocus()


###################################################
### code chunk number 37: proj2b
###################################################
my.levelplot(proj2b, x2, at = pretty(x2, 16), aspect = 1/2,
      xlab = "", ylab = "", main = "Mollweide projection")
trellis.focus("panel", 1, 1, highlight = FALSE)
print(llines(lattice.to.plot(latt, "mollweide", 30), col = 1))
trellis.unfocus()


###################################################
### code chunk number 38: interface.Rtex:591-597 (eval = FALSE)
###################################################
## stack = inla.stack(data = list(...),
##                    A = list(A1, A2, ...),
##                    effects = list(list(...),
##                                   list(...),
##                                   ...),
##                    tag = ...)


###################################################
### code chunk number 39: interface.Rtex:602-603 (eval = FALSE)
###################################################
## stack = inla.stack(stack1, stack2, ...)


###################################################
### code chunk number 40: interface.Rtex:647-651
###################################################
A = inla.spde.make.A(mesh,
                     loc = points,
                     index = rep(1:m, times = 2),
                     repl = rep(1:2, each = m) )


###################################################
### code chunk number 41: interface.Rtex:669-672
###################################################
x = as.vector(x)
covariate = rnorm(m * 2)
y = 5 + covariate * 2 + as.vector(A %*% x) + rnorm(m * 2) * 0.1


###################################################
### code chunk number 42: interface.Rtex:716-719
###################################################
mesh.index = inla.spde.make.index(name = "field",
                                  n.spde = spde$n.spde,
                                  n.repl = 2)


###################################################
### code chunk number 43: interface.Rtex:728-733
###################################################
st.est = inla.stack(data = list(y = y),
                    A = list(A, 1),
                    effects = list(c(mesh.index, list(intercept = 1)),
                                   list(cov = covariate)),
                    tag = "est")


###################################################
### code chunk number 44: interface.Rtex:753-757
###################################################
st.pred = inla.stack(data = list(y = NA),
                     A = list(1),
                     effects = list(c(mesh.index, list(intercept = 1))),
                     tag = "pred")


###################################################
### code chunk number 45: interface.Rtex:761-762
###################################################
stack = inla.stack(st.est, st.pred)


###################################################
### code chunk number 46: interface.Rtex:793-800
###################################################
formula =
  y ~ -1 + intercept + cov + f(field, model = spde, replicate = field.repl)
inla.result = inla(formula,
                   data = inla.stack.data(stack, spde = spde),
                   family = "normal",
                   control.predictor = list(A = inla.stack.A(stack),
                                            compute = TRUE))


###################################################
### code chunk number 47: postrange
###################################################
result = inla.spde2.result(inla.result, "field", spde)
plot(result[["marginals.range.nominal"]][[1]], type = "l",
     main = "Nominal range, posterior density")


###################################################
### code chunk number 48: postvar
###################################################
plot(result[["marginals.variance.nominal"]][[1]], type = "l",
     main = "Nominal variance, posterior density")


###################################################
### code chunk number 49: interface.Rtex:834-841
###################################################
index = inla.stack.index(stack, "pred")$data
linpred.mean = inla.result[["summary.linear.predictor"]]$mean
linpred.sd = inla.result[["summary.linear.predictor"]]$sd
image(proj$x, proj$y, inla.mesh.project(proj,
      linpred.mean[index[mesh.index$field.repl == 1]]))
image(proj$x, proj$y, inla.mesh.project(proj,
      linpred.sd[index[mesh.index$field.repl == 1]]))


###################################################
### code chunk number 50: interface.Rtex:848-851
###################################################
## Data -vs- linear model estimate
index.est = inla.stack.index(stack, "est")$data
plot(inla.result$summary.fitted.values$mean[index.est], y)


###################################################
### code chunk number 51: predmean
###################################################
linpred.mean = inla.result[["summary.linear.predictor"]]$mean
xplot = linpred.mean[index[mesh.index$field.repl == 1]]
my.levelplot(proj, xplot, at = pretty(xplot, 16), aspect = "iso",
      xlab = "", ylab = "", main = "", xlim = c(0, 1), ylim = c(0, 1))
trellis.focus("panel", 1, 1, highlight = FALSE)
print(lpoints(points, col = 1))
trellis.unfocus()


###################################################
### code chunk number 52: predsd
###################################################
linpred.sd = inla.result[["summary.linear.predictor"]]$sd
xplot = linpred.sd[index[mesh.index$field.repl == 1]]
my.levelplot(proj, xplot, at = pretty(xplot, 10), aspect = "iso",
      xlab = "", ylab = "", main = "", xlim = c(0, 1), ylim = c(0, 1))
trellis.focus("panel", 1, 1, highlight = FALSE)
print(lpoints(points, col = 1))
trellis.unfocus()


### R code from vignette source 'tokyo.Rtex'

###################################################
### code chunk number 1: tokyo.Rtex:1-3
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("INLA")


###################################################
### code chunk number 2: tokyo.Rtex:25-28
###################################################
data("Tokyo")
knots = seq(1, 367, length = 25)
mesh = inla.mesh.1d(knots, interval = c(1, 367), degree = 2, boundary = "cyclic")


###################################################
### code chunk number 3: tokyo.Rtex:32-35
###################################################
sigma0 = 1
kappa0 = 1e-3
tau0 = 1 / (4 * kappa0^3 * sigma0^2)^0.5


###################################################
### code chunk number 4: tokyo.Rtex:38-42
###################################################
spde = inla.spde2.matern(mesh, constr = FALSE,
                         B.tau = cbind(log(tau0), 1),
                         B.kappa = cbind(log(kappa0), 0),
                         theta.prior.prec = 1e-4)


###################################################
### code chunk number 5: tokyo.Rtex:45-58
###################################################
A = inla.spde.make.A(mesh, loc = Tokyo$time)
time.index = inla.spde.make.index("time", n.spde = spde$n.spde)
stack = inla.stack(data = list(y = Tokyo$y, link = 1, Ntrials = Tokyo$n),
                   A = list(A),
                   effects = list(time.index),
                   tag = "est")
formula = y ~ -1 + f(time, model = spde)
data = inla.stack.data(stack)
result = inla(formula, family = "binomial", data = data,
              Ntrials = data$Ntrials,
              control.predictor = list(A = inla.stack.A(stack),
                                       link = data$link,
                                       compute = TRUE))


###################################################
### code chunk number 6: tokyo
###################################################
time = 1:366
index = inla.stack.index(stack, "est")$data
plot(Tokyo$time, Tokyo$y / Tokyo$n, xlab = "Day", ylab = "Probability")
lines(time, result$summary.fitted.values$mean[index])
lines(time, result$summary.fitted.values$"0.025quant"[index], lty = 2)
lines(time, result$summary.fitted.values$"0.975quant"[index], lty = 2)


###################################################
### code chunk number 7: tokyores
###################################################
## Residuals for cumulative counts:
plot(time, cumsum(Tokyo$y / Tokyo$n) - cumsum(result$summary.fitted.values$mean[index]))


### R code from vignette source 'examples.Rtex'

###################################################
### code chunk number 1: examples.Rtex:1-5
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("INLA")
set.seed(12345L)
inla.qsample(n = 1, Matrix(1, 1, 1), seed = 12345L)


###################################################
### code chunk number 2: examples.Rtex:34-45 (eval = FALSE)
###################################################
## knots = seq(0, 100, length = 11)
## mesh1 = inla.mesh.1d(loc = knots, degree = 2, boundary = "free")
## mesh2 = inla.mesh.2d(...)
## spde = inla.spde2.matern(mesh2, alpha = 2, ...)
## index = inla.spde.make.index("space", n.spde = spde$n.spde, n.group = mesh1$m)
## formula = y ~ -1 + f(space, model = spde, group = space.group,
##                      control.group = list(model = "ar1"))
## A = inla.spde.make.A(mesh2, loc = station.loc,
##                      index = station.id, group = time,
##                      group.mesh = mesh1)
## stack = inla.stack(data = list(y = y), A = list(A), effects = list(index))


###################################################
### code chunk number 3: examples.Rtex:57-68 (eval = FALSE)
###################################################
## stack1 = inla.stack(data = list(Y = cbind(y1, NA), link = 1, N = N1), ...)
## stack2 = inla.stack(data = list(Y = cbind(NA, y2), link = 2, E = E2), ...)
## stack = inla.stack(stack1, stack2)
## data = inla.stack.data(stack)
## formula = Y ~ ...
## result = inla(formula, data = data, family = c("binomial", "poisson"),
##               Ntrials = data$N,
##               E = data$E,
##               control.predictor = list(A = inla.stack.A(stack),
##                                        link = data$link,
##                                        compute = TRUE))


