#### Test INLA and MCMC on a toy example 
library(INLA)
library(lattice)
library(fields)

#### 1 Test on the plane (example from JSS paper)
#### 1.1 Generate the True process and the data
## Simulate points -- uniform on [0,1] * [0,1]
m = 400
yloc <- matrix(runif(m*2), m, 2)
mesh <- inla.mesh.2d(loc = yloc, cutoff = 0.05, offset = c(0.1, 0.4), max.edge = c(0.05, 0.2))

## Given mean on a fine grid
xmu <- log(mesh$loc[,1]+5) - mesh$loc[,2]^2*5 + mesh$loc[,1]*mesh$loc[,2]*10
proj <- inla.mesh.projector(mesh, dims = c(100,100))
levelplot(row.values = proj$x, column.values = proj$y, x = inla.mesh.project(proj, xmu), 
          contour = TRUE, aspect = "fill", labels = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))

## simulate x ~ GP(0, Sigma(0.4, 4))
range0 <- 0.4
sigma0 <- 1
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, 1),
                          B.kappa = cbind(log(kappa0), 0, -1), 
                          theta.prior.mean = c(0,0), theta.prior.prec = c(0.1, 1))
Q <- inla.spde.precision(spde, theta = c(log(8),0))
xSim <- as.numeric(inla.qsample(n=1, Q))
## Plot the simulated process
levelplot(row.values = proj$x, column.values = proj$y, x = inla.mesh.project(proj, xSim), 
          contour = TRUE, aspect = "fill", labels = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))

## The true phyical process X = xmu + xSim
X <- xmu + xSim
## plot the true process
levelplot(row.values = proj$x, column.values = proj$y, x = inla.mesh.project(proj, X), 
          contour = TRUE, aspect = "fill", labels = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))

## Simulate the observations
## Use only a fraction of the locations
ff <- 0.7
indy <- sample(1:m, m*ff)
dataloc <- yloc[indy,]

## Find the projection matrix from mesh vertices to the locations
A <- inla.spde.make.A(mesh, loc = dataloc)

## Generate the data y = Ax +  e
errors <- rep(0, m*ff)
errors[which(dataloc[,1] <=0.25)] <- 0.02
errors[which(dataloc[,1] > 0.25 & dataloc[,1] <= 0.5)] <- 0.05
errors[which(dataloc[,1] > 0.5 & dataloc[,1] <= 0.75)] <- 0.1
errors[which(dataloc[,1] > 0.75)] <- 0.5
y <- as.vector(A %*% X) + rnorm(m*ff)*errors

#### 1.2 Now start INLA analysis
ydata <- y - as.vector(A %*% xmu)

st.est <- inla.stack(data = list(y=ydata), A = list(A),
                     effects = list(GIA = 1:spde$n.spde), tag = "est")

## Predict at the y location, mesh grid and predict error location
xg <- seq(0, 1, length.out = 10)
yg <- seq(0, 1, length.out = 10)
xyg <- as.matrix(expand.grid(x = xg, y = yg))
A_pred <- inla.spde.make.A(mesh = mesh, loc = xyg)

st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(GIA=1:spde$n.spde), tag = "pred")

stGIA <- inla.stack(st.est, st.pred)


hyper <- list(prec = list(fixed = TRUE, initial =0))
formula = y ~ -1 + f(GIA, model = spde)
prec_scale <- c(1/errors^2, rep(1, 100))
res_inla <- inla(formula, data = inla.stack.data(stGIA, spde = spde), family = "gaussian",
                 scale = prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stGIA), compute =TRUE))
summary(res_inla)

## Plot the posteriors of the parameters
result <- inla.spde2.result(res_inla, "GIA", spde)
par(mfrow= c(1,2))
plot(result[["marginals.range.nominal"]][[1]], type = "l",
     main = "Nominal range, posterior density")
plot(result[["marginals.variance.nominal"]][[1]], type = "l",
     main = "Nominal variance, posterior density")

## Plot the predicted GIA field mean and variance
pidx <- inla.stack.index(stGIA, tag = "pred")
GIA_mpost <- res_inla$summary.random$GIA$mean + xmu
GIA_spost <- res_inla$summary.random$GIA$sd

xyg_mpost <- res_inla$summary.linear.predictor$mean[pidx$data]
xyg_spost <- res_inla$summary.linear.predictor$sd[pidx$data]

y_mpost <- res_inla$summary.linear.predictor$mean[1:(m*ff)] + A%*%xSim
y_spost <- res_inla$summary.linear.predictor$sd[1:(m*ff)]


## Plot the variance
par(mfrow = c(3,1))
theta_mode <- result$summary.theta$mode
Q2 <- inla.spde2.precision(spde, theta = theta_mode)
Q1 <- inla.spde.precision(spde, theta = c(0,0))
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(sqrt(1/diag(Q1)))), col = topo.colors(7),
           xlim = c(0,1), ylim = c(0,1), breaks = c(0, 0.02, 0.05, 0.1,0.5, 2, 4, 8),
           xlab = "Longitude", ylab = "Latitude", main = "Matern prior Error field")
points(dataloc[,1], dataloc[,2], cex = y_spost*2, pch = 1)

image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(sqrt(1/diag(Q2)))), col = topo.colors(7),
           xlim = c(0,1), ylim = c(0,1), breaks = c(0, 0.02, 0.05, 0.1,0.5, 2, 4, 8),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field")
points(dataloc[,1], dataloc[,2], cex = y_spost*2, pch = 1)


## The standard error
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_spost)), col = topo.colors(7),
           xlim = c(0,1), ylim = c(0,1), breaks = c(0, 0.02, 0.05, 0.1,0.5, 2, 4, 8),
           xlab = "Longitude", ylab = "Latitude", main = "predicted error field")
points(dataloc[,1], dataloc[,2], cex = y_spost*2, pch = 1)
points(dataloc[,1], dataloc[,2], cex = errors*2, pch = 1, col = 2)

## Plot the error on mean
par(mfrow = c(3,1))
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(X)), col = topo.colors(40),
           breaks = seq(-10, 10, 0.5),
           xlab = "Longitude", ylab = "Latitude", main = "True Process",
           xlim = c(0, 1), ylim = c(0,1))
image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(xmu)), col = topo.colors(40),
           xlab = "Longitude", ylab = "Latitude", main = "The prior mean",
           breaks = seq(-10, 10, 0.5),
           xlim = c(0, 1), ylim = c(0,1))

image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIA_mpost)), col = topo.colors(40),
           breaks = seq(-10, 10, 0.5),
           xlab = "Longitude", ylab = "Latitude", main = "Matern posterior Error field",
           xlim = c(0, 1), ylim = c(0,1))
points(xyg[,1], xyg[,2], cex = xyg_spost)
points(dataloc[,1], dataloc[,2], cex = y_spost, pch = 1, col = 2)


## Add some signals to the data
## 15% points nearer southwest are lower
n1 = round(m*ff*0.1)
id1 <- sample(which(yyloc[,1] <= 0.5 & yyloc[,2] <= 0.5), n1)
y[id1] <- y[id1] - abs(rnorm(n1, sd = 3))
levelplot(row.values = proj$x, column.values = proj$y, x = inla.mesh.project(proj, xSim), 
          panel=function(...){
            panel.levelplot(...)
            panel.points(yyloc[id1,], col = "red", pch = 19)
          },
          contour = TRUE, aspect = "fill", labels = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))

## 5% points nearer north are higher
n2 <- round(m*ff*0.10)
id2 <- sample(which(yyloc[,2] > 0.7), n2)
y[id2] <- y[id2] + abs(rnorm(n2, sd = 2))
levelplot(row.values = proj$x, column.values = proj$y, x = inla.mesh.project(proj, xSim), 
          panel=function(...){
            panel.levelplot(...)
            panel.points(yyloc[id1,], col = "red", pch = 19)
            panel.points(yyloc[id2,], col = "blue", pch = 19)
          },
          contour = TRUE, aspect = "fill", labels = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))


