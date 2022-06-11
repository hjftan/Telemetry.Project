library(sp)
library(raster)
library(ggplot2)
library(tidyverse)
library(sf)
library(quantmod)
library(rgeos)
library(INLA)
library(inlabru)
library(prioritizr)
##

load("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Telemetry.create.extract.enviro.data.RData")

################################################################################################################################################
#this script 
#- uses output of script Telemetry.create.extract.enviro.data.R
#- fits a bayesian spatiotemporal mixed effects model to movement data from Caspian seals (Pusa caspica)

#- the seal data was collected using Argos satellite telemetry tags and processed
#using a state space model
#- From the processed data I calculate the movement persistence metric https://github.com/ianjonsen/mpmm 
#and associate foraging-like behaviors with Sea surface temperature (SST), proximity to river inlets, 
#and proximity to the 50m isobath bathymetric contour
#a fixed effect for SST is defined because the relationship between SST and foraging behavior is 
#assumed to be constant throughout the year
#fixed effects are also defined for proximity to river inlets and proximity to the 50m isobath
#temporal deviations are also included for proximity to river inlets and proximity to the 50m isobath
#because foraging activity near these features is likely to be driven by prey presence
#and the fish prey of Caspian seals are believed to exhibit spatially and temporally structured
#lifecycles
#- individual based random temporal effects and global spatiotemporal effects are also included

################################################################################################################################################

caspian.shape <- spTransform(caspian.shape, CRS(proj4string(all.data)))

##################################################################
#creates a 2 dimensional mesh to define a Stochastic Partial Differential Equation (SPDE) based
#spatial effect

mesh.c <- 50
mesh.dis <- mesh.c * 2

mesh_bnd.2 <- spsample(caspian.shape, n = 10000, type = "regular")
mesh_bnd.1 <- INLA::inla.nonconvex.hull(mesh_bnd.2, convex = mesh.c, concave = mesh.c * 10)
mesh_bnd.2 <- INLA::inla.nonconvex.hull(mesh_bnd.2, convex = mesh.c)#, concave = mesh.c)

mesh <- INLA::inla.mesh.2d( # SB
  boundary = mesh_bnd.1,
  max.edge = mesh.c,
  #min.angle = c(25,19),
  #cutoff = 100,
  crs = CRS(proj4string(all.data))
)

my.pcmatern <- inla.spde2.pcmatern(mesh,
                                   prior.sigma = c(1, 0.01), # P(??>??0)=ps.
                                   prior.range = c(100, 0.01), constr = FALSE) # P(r<r0)=pr. # c(150, 0.50))

#####################################################################################################################
#creates a 1 dimensional meshs to define a SPDE based temporal effect

mesh.seq.temp <- seq(from = min(all.data$idate), to = max(all.data$idate), length.out = 12)
mesh.seq.cov <- seq(from = min(all.data$idate), to = max(all.data$idate), length.out = 12)

mesh1D.idate.free <- inla.mesh.1d(loc = mesh.seq.temp, boundary = "free", degree = 2)
mesh1D.idate.cyc <- inla.mesh.1d(loc = mesh.seq.cov, boundary = "cyclic", degree = 2)
#free for idate as the mean of the effect is estimated from each individual
#and most tags do not last close to a full calender year
#cyclic for temporal as the mean of this effect is estimated from the full dataset
#and we have good coverage across a ful calender year

mat1D.idate.temp <- inla.spde2.pcmatern(mesh1D.idate.free,
                                        prior.range = c(10, 0.01), prior.sigma = c(1, 0.01),
                                        constr = FALSE)
mat1D.idate.covs <- inla.spde2.pcmatern(mesh1D.idate.cyc, 
                                        prior.range = c(50, 0.01), prior.sigma = c(1, 0.01), 
                                        constr = TRUE)
#units of range on these effects are days. define the prior for the average range on the temporal effect
#as 10 as individual variation in behaviours can be variable
#a higher value for range is used for the temporal covariate effect as this is likely going to be
#driven by broader ecosystem changes and not individual variations in behaviour

#integrate to zero constraints are included on the covairate effects so the full effect of each
#convariate is equal to the fixed effect + temporal effect
#they are not included on the indiviudal based temporal effect so this effect represents the individual
#based mean for the temporal trend in the behavioural metric

#######################################################################

all.data$index.idate <- inla.group(all.data$idate, n = 12, method = "cut", idx.only = TRUE)
mesh1D.index.idate <- inla.mesh.1d(loc = 1:max(all.data$index.idate), boundary = "free", degree = 2)

###############################################################################################################

cmp.time <- ~
  -1 +
  idv.smooth(main = idate, model = mat1D.idate.temp,
             group = index.id, ngroup = max(all.data$index.id),
             control.group = list(model = "iid")) +
  #individual based temporal effect with a integrate to zero constaint
  riv.cov(cent.sq.riv, model = "linear", scale.model = FALSE, mean.linear = 0, prec.linear = 1) +
  con.cov(cent.sq.dis.con, model = "linear", scale.model = FALSE, mean.linear = 0, prec.linear = 1) +
  sst.cov(cent.sst, model = "linear", scale.model = FALSE, mean.linear = 0, prec.linear = 1) +
  #fixed effects for environmental covariates
  riv.cov.time(idate, cent.sq.riv, model = mat1D.idate.covs) +
  con.cov.time(idate, cent.sq.dis.con, model = mat1D.idate.covs) +
  #temporal effects (deviations from the fixed effect) are included for rivers and contour
  spat.field(main = coordinates, model = my.pcmatern,
             group = index.idate, group_mapper = bru_mapper(mesh1D.index.idate, indexed=TRUE),
             control.group = list(model="ar1"))
#spatialtemporal SPDE effect to account for large amounts of unexplained variation in siteuse/foraging activity


form.time <- g.auc.stan.inv ~
  idv.smooth +
  
  riv.cov +
  con.cov +
  sst.cov +

  riv.cov.time +
  con.cov.time +

  spat.field

model.time <- bru(components = cmp.time,
                  formula = form.time,
                  data = all.data,
                  family = "beta",
                  options = list(bru_verbose = 3,
                                 inla.mode="experimental",
                                 control.inla = list(h = 0.005,
                                                     int.strategy = "eb"),
                                 bru_max_iter = 5))


#just to check if having an explicit inter does anything in our case
#save.image("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Telemetry.fit.inla.model.RData")

