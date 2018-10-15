#####################################################
########  Plot and analysis the MCMC results ########
#####################################################

library(coda)
library(lattice)
library(fields)

theta_1 <- mcmc.list(lapply(res, function(x) mcmc(x$samples$theta1, start = burnin+1,
                                                  end = burnin+numsamples*thinning, thin = thinning)))
theta_2 <- mcmc.list(lapply(res, function(x) mcmc(x$samples$theta2, start = burnin+1,
                                                  end = burnin+numsamples*thinning, thin = thinning)))
xGIA <- list()
for (i in c(1, 50)){
  xGIA[[i]] <- mcmc.list(lapply(res, function(x) mcmc(x$samples$xGIA[,i], start = burnin+1,
                                                      end = burnin+numsamples*thinning, thin = thinning)))
}

#### plot Sample analysis
## traceplot
pdf(file = paste0(wkdir, exname, "_MCMCanalysis.pdf"), width = 8, height = 6)
par(mfrow = c(2,2))
traceplot(theta_1, main = expression(theta[1]))
traceplot(theta_2, main = expression(theta[2]))

traceplot(xGIA[[1]], main = expression(x[1]))
traceplot(xGIA[[50]], main = expression(x[50]))

print(acfplot(theta_1, main = expression(theta[1])))
print(acfplot(theta_2, main = expression(theta[2])))

print(acfplot(xGIA[[1]], main = expression(x[1])))
print(acfplot(xGIA[[50]], main = expression(x[50])))


print(densityplot(theta_1, main = expression(theta[1])))
print(densityplot(theta_2, main = expression(theta[2])))

print(densityplot(xGIA[[1]], main = expression(x[1])))
print(densityplot(xGIA[[50]], main = expression(x[50])))

dev.off()

#### The GIA field
## Posterior mean
GIApost_mean <- lapply(res, function(x) colMeans(x$samples$xGIA))
GIApost_sd <- lapply(res, function(x) apply( x$samples$xGIA, 2, sd))
GIA_ice6g2 <- cbind(GIA_ice6g, GIApost_mean, GIApost_sd)

proj <- inla.mesh.projector(Mesh_GIA, projection = "longlat", dims = c(361,181))

for (i in (1:n.chains)){
  pdf(file = paste0(wkdir, exname, "_GIAfield", i, ".pdf"), width = 12, height = 8)
  par(mfrow = c(2,1))
  image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIApost_mean[[i]])), col = topo.colors(40),
             xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- mean")
  points(GPSX, GPSY, pch = 20)

  image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(GIApost_sd[[i]])), col = topo.colors(40),
             xlab = "Longitude", ylab = "Latitude", main = "Posterier marginals -- sd")
  points(GPSX, GPSY, pch = 20)
  dev.off()
}

pdf(file = paste0(wkdir, exname, "_hyperpars.pdf"), width = 8, height = 4)
for (i in (1:n.chains)){
  range_samp <- exp(lmu_r + theta_2[[i]])
  dens_rho <- density(range_samp)
  rho_mode <- dens_rho$x[dens_rho$y==max(dens_rho$y)]

  sigma_samp <- exp(lmu_s + theta_1[[i]])
  dens_s <- density(sigma_samp)
  sigma_mode <- dens_s$x[dens_s$y==max(dens_s$y)]

  par(mfrow=c(1,2))
  plot(dens_rho, main = bquote(bold(rho("mode") == .(round(rho_mode, 4)))))
  plot(dens_s,  main = bquote(bold(sigma("mode") == .(round(sigma_mode, 4)))))

}
dev.off()


GIA_samps <- do.call(rbind, lapply(res, function(x) x$samples$xGIA))
GIA_mean <- colMeans(GIA_samps)
GIApost_cov <- cov(GIA_samps)

GIA_mpred <- as.numeric(CMat2 %*% GIA_mean)
CMat2C <- CMat2 %*% GIApost_cov
GIA_spred <- rep(0, nrow(CMat2))

ids <- c(seq(1, 64800, 6480), 64800)
GIA_spred <- mclapply.hack(1:10, function(i) sapply(ids[i]:ids[i+1], function(x) CMat2[x,] %*% CMat2C[x,]))
GIA_spreds <- c(sapply(1:10, function(x) GIA_spred[[x]][1:6480]))
GIA_ice6g2 <- cbind(GIA_ice6g, GIA_mpost = GIA_mpred, GIA_spost = GIA_spreds)
write.table(GIA_ice6g2, file = paste0(wkdir, exname, "_GIApred.txt"), row.names = FALSE, eol = "\r\n")
