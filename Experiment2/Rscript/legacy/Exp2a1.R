
## ----loadlibs, message = FALSE, warning = FALSE--------------------------
library(rgdal); library(sp);library(GEOmap)
library(INLA)
source("~/glbm/BHM_sphere/functions.R")
source("~/glbm/BHM_sphere/partition_fun.R")

## ----alt_data------------------------------------------------------------
library(ncdf4)
alt_nc <- nc_open("/./projects/GlobalMass/WP3-Ocean/BHMinputs/SSH/trend_SSH_CCI_200501_201512.nc")
print(alt_nc)
lon <- ncvar_get(alt_nc, "lon")
lat <- ncvar_get(alt_nc, "lat")
trend_ssh <- ncvar_get(alt_nc, "trend_ssh") #note that there re NAs for land datat
err_ssh <- ncvar_get(alt_nc, "err")

alt_data <- data.frame(trend_ssh = as.numeric(trend_ssh), err_ssh = as.numeric(err_ssh))
alt_data$lon <- rep(lon, 180) 
alt_data$lat <- rep(lat, each = 360)
alt_data2 <- na.omit(alt_data)
## Find the xyz coords of the altimetry data
alt_loc <- do.call(cbind, Lll2xyz(lon = alt_data2$lon, lat = alt_data2$lat))


## ----prior---------------------------------------------------------------
## Priors mean and variance for the parameters: rho and sigma
mu_r <- 3000/6371
v_r <- 1
mu_s <- 5
v_s <- 10^2

## ----ssh_mesh------------------------------------------------------------
## Genereate Fibonacci points on the sphere
fibo_points <- fiboSphere(N = 12960, L0 = TRUE)
fibo_points_xyz <- do.call(cbind, Lll2xyz(lat = fibo_points[,2], lon = fibo_points[,1]))
mesh0 <- inla.mesh.2d(loc = fibo_points_xyz, cutoff = 0.01, max.edge = 1)
## Make this "smoother"
mesh0 <- inla.mesh.2d(loc = mesh0$loc, cutoff = 0.01, max.edge = 1)

## ----ssh_mesh2-----------------------------------------------------------
## Load the Ocean polygon
Ocean <- readOGR(dsn = "/./projects/GlobalMass/WP1-BHM/maps/ne_110m_ocean", layer = "ne_110m_ocean")

## Remove mesh not in the ocean
mesh0 <- dt.mesh.addon.posTri(mesh = mesh0, globe = TRUE) 
Tlonlat <- Lxyz2ll(list(x = mesh0$posTri[,1], y = mesh0$posTri[,2], z = mesh0$posTri[,3]))
mesh0$Trill <- cbind(lon = Tlonlat$lon, lat =Tlonlat$lat)
TinOcean <- unlist(over(Ocean, SpatialPoints(coords=mesh0$Trill, proj4string = CRS(proj4string(Ocean))), returnList=T))
TAll <- 1:mesh0$t
ToutOcean <- TAll[-TinOcean]
Omega = dt.Omega(list(TinOcean, 1:mesh0$t), mesh0)

mesh_ssh <- mesh.sub(mesh0, Omega, 1)
summary(mesh_ssh)

## ----spde----------------------------------------------------------------
## Transform the parameters for the SPDE_GMRF approximation
trho <- Tlognorm(mu_r, v_r)
tsigma <- Tlognorm(mu_s, v_s)

## Build the SPDE model with the prior
lsigma0 <- tsigma[1]
theta1_s <- tsigma[2]
lrho0 <- trho[1]
theta2_s <- trho[2]
lkappa0 <- log(8)/2 - lrho0
ltau0 <- 0.5*log(1/(4*pi)) - lsigma0 - lkappa0

SSH_spde <- inla.spde2.matern(mesh_ssh, B.tau = matrix(c(ltau0, -1, 1),1,3), B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                              theta.prior.mean = c(0,0), theta.prior.prec = c(sqrt(1/theta1_s), sqrt(1/theta2_s)))


## ----link_data-----------------------------------------------------------
## Link the process to observations and predictions
A_data <- inla.spde.make.A(mesh = mesh_ssh, loc = alt_loc)
A_pred <- inla.spde.make.A(mesh = mesh_ssh, loc = rbind(mesh_ssh$loc))

## Create the estimation and prediction stack
st.est <- inla.stack(data = list(y=alt_data2$trend_ssh), A = list(A_data),
                     effects = list(SSH = 1:SSH_spde$n.spde), tag = "est")
st.pred <- inla.stack(data = list(y=NA), A = list(A_pred),
                      effects = list(SSH=1:SSH_spde$n.spde), tag = "pred")
stSSH <- inla.stack(st.est, st.pred)

## ----inla_run, include = TRUE, eval = FALSE------------------------------
## Fix altimetry errors as they are known
hyper <- list(prec = list(fixed = TRUE, initial = 0))
prec_scale <- c(1/alt_data2$err_ssh^2, rep(1, nrow(A_pred)))

## The formular for modelling the SSH mean
formula = y ~ -1 +  f(SSH, model = SSH_spde)

## Run INLA
res_inla <- inla(formula, data = inla.stack.data(stSSH, spde = SSH_spde), family = "gaussian",
                 scale =prec_scale, control.family = list(hyper = hyper),
                 control.predictor=list(A=inla.stack.A(stSSH), compute =TRUE))

## ----inla_res, include = TRUE, eval = FALSE------------------------------
INLA_pred <- res_inla$summary.linear.predictor
## Extract and project predictions
pred_idx <- inla.stack.index(stSSH, tag = "pred")$data

## SSH
SSH_m <- INLA_pred$mean[pred_idx]
SSH_u <- INLA_pred$sd[pred_idx]
proj <- inla.mesh.projector(mesh_ssh, projection = "longlat", dims = c(360,180), xlim = c(0,359), ylim = c(-89.5, 89.5))
SSH_grid <- expand.grid(proj$x, proj$y)
SSH_pred <- data.frame(lon = SSH_grid[,1], lat = SSH_grid[,2],
                       mean = as.vector(inla.mesh.project(proj, as.vector(SSH_m))),
                       u = as.vector(inla.mesh.project(proj, as.vector(SSH_u))))

res_SSH <- list(res_inla = res_inla, spde = SSH_spde, st = stSSH,
            mesh = mesh_ssh,  SSH_pred = SSH_pred)

save(res_SSH, file ="/./projects/GlobalMass/WP1-BHM/Experiment2a/exp2a_ssh.RData")
