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

load("Telemetry.fit.inla.model.RData")

################################################################################################################################################
#this script 
#- uses output of script Telemetry.fit.inla.model.R
#- predicts foraging activity across space and time
#- generates marine protected area plans (MPA) using a prioritization algorithm
#- creates some summary outputs of vessel activity within MPAS

########################################################################################################################################
#read in vessel rasters

read.ras <- function(){
  for ( i in 1:12){
    #ras <- raster(paste0("vessel.rasters/season.num.ves.", i, ".grd"))
    #assign(paste0("season.num.ves.", i), ras, envir = globalenv())
    ras <- raster(paste0("vessel.rasters/season.num.ves.", i, "new.grd"))
    assign(paste0("season.num.ves.", i), ras, envir = globalenv())
  }
} 
read.ras() 

vessel.data.all <- stack(season.num.ves.1, season.num.ves.2, season.num.ves.3, season.num.ves.4,
                         season.num.ves.5, season.num.ves.6, season.num.ves.7, season.num.ves.8,
                         season.num.ves.9, season.num.ves.10, season.num.ves.11, season.num.ves.12)
# rasters are number of vessels in each season per cell

vessel.data.all <- vessel.data.all/
  (mesh.seq.temp[2] - mesh.seq.temp[1]) 
# num vessels per day

vessel.data.all <- vessel.data.all/res(vessel.data.all)[1] 
# num vessels per day per km

############################################################################################################################################
#do spatial predictions

date.titles <- as.list(1:max(fit.mpm$index.idate))
date.titles <- lapply(date.titles, function(x){
  dates <- subset(fit.mpm, index.idate == x)$date
  
  min.d <- as.Date(paste0(1966, substr(dates, 5, 10)), tryFormats = c("%Y-%m-%d"))
  max.d <- as.Date(paste0(1967, substr(dates, 5, 10)), tryFormats = c("%Y-%m-%d"))
  
  if((sum(grepl("12", substr(max.d, 6, 7)) != 0)) & (sum(grepl("01", substr(max.d, 6, 7)) != 0))){
    min.d <- min.d[-which(as.numeric(substr(min.d, 6, 7)) < 6)]
    max.d <- max.d[-which(as.numeric(substr(max.d, 6, 7)) > 6)]
  }
  
  range <- range(min.d, max.d) %>%
    substr(., 6, 10) %>%
    paste(., sep=" ", collapse=" to ") 
  
  return(range)
})
date.titles[[6]] <- "02-12 to 04-16" # not data between 04-13 and 04-16 to manually change thuis


fit.mpm$iseason.ras <- inla.group(fit.mpm$idate, n = 12, idx.only = TRUE) # same as iseason# see hpc

#create dataframe with indexs for predctions which can be matched to vessel data
match.index.idate <- 
  unique(fit.mpm@data[,c("iseason.ras", "index.idate", "idate", "date")]) %>%
  split(., .$iseason.ras) %>%
  lapply(., function(x){
    return(data.frame(iseason.ras = unique(x$iseason.ras),
                      index.idate = mean(x$index.idate), 
                      idate = mean(x$idate),
                      minidate = min(x$idate),
                      maxidate = max(x$idate)))}) %>%
  do.call("rbind", .)

spat.pred <- crop(as(vessel.data.all[[1]], "SpatialPixels"),caspian.shape)
spat.pred <- SpatialPixelsDataFrame(spat.pred, data.frame(dumb = rep(1, length(spat.pred))))

spat.pred <- cprod(spat.pred, match.index.idate)
rm(match.index.idate)
#######################################################################################################################################
#extact environmnetal data for the space where we want t o generate spatial predictions and scale to model
spat.pred$Dis.Any <- raster::extract(distance.any.ras, spTransform(as(spat.pred, "SpatialPointsDataFrame"),
                                                                   CRS(proj4string(distance.any.ras))))
spat.pred$dis.con <- raster::extract(contour, spTransform(as(spat.pred, "SpatialPointsDataFrame"),
                                                          CRS(proj4string(contour))))
spat.pred$dis.haul.any <- raster::extract(dis.haulouts, spTransform(as(spat.pred, "SpatialPointsDataFrame"),
                                                                    CRS(proj4string(dis.haulouts))))
spat.pred$bath <- raster::extract(Bathymetry, spTransform(as(spat.pred, "SpatialPointsDataFrame"),
                                                          CRS(proj4string(Bathymetry))))

unique(fit.mpm$index.idate)
unique(fit.mpm$date)

data.ind.df <- unique(fit.mpm@data[,c("index.idate", "date", "idate", "iseason.ras")])
data.ind.df$date <- gsub("-", "", data.ind.df$date)

SST <- c(list.files("SST.MUR.NEW/2008", full.names = T, recursive = T),
         list.files("SST.MUR.NEW/2009", full.names = T, recursive = T),
         list.files("SST.MUR.NEW/2010", full.names = T, recursive = T),
         list.files("SST.MUR.NEW/2011", full.names = T, recursive = T),
         list.files("SST.MUR.NEW/2012", full.names = T, recursive = T),
         #list.files("SST.MUR.NEW/2013", full.names = T, recursive = T),
         list.files("SST.MUR.NEW/2016", full.names = T, recursive = T),
         list.files("SST.MUR.NEW/2017", full.names = T, recursive = T))#,
#list.files("SST.MUR.NEW/2018", full.names = T, recursive = T))
#sst data from the MUR project

SST.list <- list()
for ( i in 1:max(fit.mpm$iseason.ras)){
  my.sst <- stack(unique(grep(paste(subset(data.ind.df, iseason.ras == i)$date, collapse="|"),
                              SST, value=TRUE)))
  my.sst <- mean(my.sst, na.rm = T)
  SST.list[[i]] <- my.sst
}
rm(my.sst)
SST.list

extract.sst.data <- function(data){
  
  data <- split(data, data$index.idate)
  data <- lapply(data, function(x){
    
    my.sst <- SST.list[[unique(x$index.idate)]]
    
    my.crs <- CRS(proj4string(x))
    x.p <- as(x, "SpatialPointsDataFrame")
    x.p <- spTransform(x.p, CRS(proj4string(my.sst)))
    
    x$SST <- raster::extract(my.sst, x.p)
    
    return(x)
  }
  )
  data <- do.call("rbind", data)
  return(data)
}
spat.pred <- extract.sst.data(spat.pred)
spat.pred <- subset(spat.pred, !is.na(SST)) 

spat.pred$cent.sst <- (spat.pred$SST - mean(c(fit.mpm$SST))) / sd(c(fit.mpm$SST))
spat.pred$sq.riv <- sqrt(spat.pred$Dis.Any)
spat.pred$cent.sq.riv <- (spat.pred$sq.riv - mean(c(fit.mpm$sq.riv))) / sd(c(fit.mpm$sq.riv))

spat.pred$sq.dis.con <- sqrt(spat.pred$dis.con)
spat.pred$cent.sq.dis.con <- (spat.pred$sq.dis.con - mean(c(fit.mpm$sq.dis.con))) / sd(c(fit.mpm$sq.dis.con))

spat.pred$sq.dis.haul.any <- sqrt(spat.pred$dis.haul.any)
spat.pred$cent.sq.dis.haul.any <- (spat.pred$sq.dis.haul.any - mean(c(fit.mpm$sq.dis.haul.any))) / sd(c(fit.mpm$sq.dis.haul.any))

spat.pred$bath.abs <- abs(spat.pred$bath)
spat.pred$cent.bath.abs <- (spat.pred$bath.abs - mean(c(fit.mpm$bath.abs))) / sd(c(fit.mpm$bath.abs))

spat.pred <- predict(model.time, spat.pred, ~ 
                       exp(riv.cov +
                             con.cov +
                             sst.cov +
                             riv.cov.time +
                             con.cov.time + 
                             spat.field
                       ) /
                       (1 + exp(riv.cov +
                                  con.cov +
                                  sst.cov +
                                  riv.cov.time +
                                  con.cov.time + 
                                  spat.field
                       )),
                     include = c("riv.cov",
                                 "con.cov",
                                 "sst.cov",
                                 "riv.cov.time",
                                 "con.cov.time", 
                                 "spat.field"
                     ))


################################################################################################################################
#creates protected area plans

resolution.fact <- 1 # # to speedup computation make more than 1

create.priority.zones <- function(spatial.pred, vessel.data, max.clustering, cover, resolution.fact){
  
  all.priority.zones <- list()
  my.priority.zones.1 <- list()
  my.priority.zones.2 <- list()
  my.priority.zones.free.1 <- list()
  my.priority.zones.free.2 <- list()
  
  max.avg.dis <- 
    as.list(1:max(fit.mpm$iseason.ras)) %>%
    lapply(., function(x){
      x <- subset(fit.mpm, iseason.ras == x)
      my.dis <- gDistance(x, byid = TRUE)
      x <- mean(my.dis) + (sd(my.dis))
      return(x)
    }) %>%
    do.call("c", .) %>%
    max(.)
  
  max.ars.sum <- 
    as.list(1:max(fit.mpm$iseason.ras)) %>%
    lapply(., function(x){
      x <- subset(fit.mpm, iseason.ras == x)
      x <- sum(x$g.auc.stan.inv)
      return(x)
    }) %>%
    do.call("c", .) %>%
    max(.)
  
  any.ones <- 
    split(fit.mpm, fit.mpm$index.id) %>%
    lapply(., function(y){return(nrow(y))}) %>%
    do.call("c", .)
  if(1 %in% any.ones){
    any.tracks <- 
      split(fit.mpm, fit.mpm$index.id) %>%
      .[-which(any.ones == 1)]
  } else{
    any.tracks <- split(fit.mpm, fit.mpm$index.id)
  }
  any.tracks <- #find cells with no seal tracks over them in any season
    any.tracks %>%
    lapply(., function(x){ # make points to lines
      my.proj <- proj4string(fit.mpm)
      x <- spTransform(x, CRS(my.proj))
      x <- x[order(x$idate, decreasing = FALSE),]
      row.names(x) <- 1:nrow(x)
      x$x.coord <- x@coords[,1]
      x$y.coord <- x@coords[,2]
      x$to.x.coord <- c(x$x.coord[2:nrow(x)], 0) 
      x$to.y.coord <- c(x$y.coord[2:nrow(x)], 0) 
      
      rows <- x@data
      rows <- rows[(-nrow(rows)),c("x.coord", "y.coord", "to.x.coord", "to.y.coord")]
      
      rows <- split(rows, seq(nrow(rows)))
      
      lines <- lapply(rows, function(row) {
        lmat <- matrix(unlist(row), ncol = 2, byrow = TRUE)
        return(st_linestring(lmat))
      })
      lines <- st_sfc(lines)
      lines <- as(lines, "Spatial")
      
      proj4string(lines) <- CRS(my.proj)
      
      return(lines)
    }) %>% # make points to lines
    do.call("rbind", .) %>%
    rasterize(x = .,
              y = raster(spatial.pred),
              fun = function(y){
                return(length(y[!is.na(y)]))
              })
  
  any.tracks <- any.tracks>0
  any.tracks <- abs(any.tracks - 1)
  any.tracks <- aggregate(any.tracks, fact = resolution.fact)
  
  binary.ves <- mean(vessel.data)
  binary.ves <- binary.ves > (mean(binary.ves[]))
  binary.ves <- aggregate(binary.ves, fact = resolution.fact)
  
  for(i in 1:max(spat.pred$iseason.ras)){
    
    #using a fixed area now. fixed is more informative for the comaprison with free
    #my.ars.sum <- sum(subset(fit.mpm, iseason.ras == i)$g.auc.stan.inv)
    #sim.targ <- base.cover * (my.ars.sum / max.ars.sum)
    
    my.dis <- gDistance(subset(fit.mpm, iseason.ras == i), byid = TRUE)
    my.avg.dis <- mean(my.dis) + (sd(my.dis))
    my.boundary.penalty <- max.clustering * (my.avg.dis / max.avg.dis)
    
    any.ones <- 
      subset(fit.mpm, iseason.ras == i) %>%
      split(., .$index.id) %>%
      lapply(., function(y){return(nrow(y))}) %>%
      do.call("c", .)
    
    if(1 %in% any.ones){
      my.tracks <- 
        subset(fit.mpm, iseason.ras == i) %>%
        split(., .$index.id) %>%
        .[-which(any.ones == 1)]
    } else{
      my.tracks <- 
        subset(fit.mpm, iseason.ras == i) %>%
        split(., .$index.id)
    }
    my.tracks <- #find cells with no track over them. from current season
      my.tracks %>%
      lapply(., function(x){ # make points to lines
        
        my.proj <- proj4string(fit.mpm)
        x <- spTransform(x, CRS(my.proj))
        
        x <- x[order(x$idate, decreasing = FALSE),]
        row.names(x) <- 1:nrow(x)
        x$x.coord <- x@coords[,1]
        x$y.coord <- x@coords[,2]
        x$to.x.coord <- c(x$x.coord[2:nrow(x)], 0) 
        x$to.y.coord <- c(x$y.coord[2:nrow(x)], 0) 
        
        rows <- x@data
        rows <- rows[(-nrow(rows)),c("x.coord", "y.coord", "to.x.coord", "to.y.coord")]
        
        rows <- split(rows, seq(nrow(rows)))
        lines <- lapply(rows, function(row) {
          lmat <- matrix(unlist(row), ncol = 2, byrow = TRUE)
          return(st_linestring(lmat))
        })
        lines <- st_sfc(lines)
        lines <- as(lines, "Spatial")
        
        proj4string(lines) <- CRS(my.proj)
        
        return(lines)
      }) %>% # make points to lines
      do.call("rbind", .) %>%
      rasterize(x = .,
                y = raster(spatial.pred),
                fun = function(y){
                  return(length(y[!is.na(y)]))
                })
    
    my.tracks <- my.tracks<=0
    
    my.tracks <- aggregate(my.tracks, fact = resolution.fact)
    my.spat.pred <- aggregate(raster(subset(spatial.pred, iseason.ras == i)[,"median"]), fact = resolution.fact) 
    my.ves.ras <- aggregate(vessel.data[[i]], fact = resolution.fact)
    
    my.tracks[which(is.na(my.spat.pred)[]==TRUE)] <- NA
    binary.ves[which(is.na(my.spat.pred)[]==TRUE)] <- NA
    my.spat.pred[which(is.na(my.spat.pred)[]==TRUE)] <- NA
    my.ves.ras[which(is.na(my.spat.pred)[]==TRUE)] <- NA
    any.tracks[which(is.na(my.spat.pred)[]==TRUE)] <- NA
    
    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####
    
    my.ves.ras <- my.ves.ras + 1 # so no zero costs
    my.free.ves <- my.ves.ras
    my.free.ves[!is.na(my.free.ves)] <- mean(my.ves.ras[], na.rm = T) #need this free bit seperate so can  still include the real one in final areas for how much plots
    
    for(covi in 1:length(cover)){
      priorit.solve <- 
        problem(my.ves.ras, my.spat.pred) %>%
        add_min_set_objective() %>% 
        add_relative_targets(targets = cover[covi]) %>% 
        
        add_linear_penalties(penalty = 5, data = my.tracks) %>% #prefer cells with real trajectory data in current season
        add_linear_penalties(penalty = 2.5, data = any.tracks) %>% #prefer cells cells with real trajectory data in any season
        add_linear_penalties(penalty = 5, data = binary.ves) %>% #prefer cells which are not regularly used routes across any season
        add_boundary_penalties(penalty = my.boundary.penalty, edge_factor = 0) %>% 
        add_gurobi_solver(time_limit = 1000, threads = 8, gap = 0.1) %>% 
        solve()
      
      priorit.solve.free <- 
        problem(my.free.ves, my.spat.pred) %>%
        add_min_set_objective() %>% 
        add_relative_targets(targets = cover[covi]) %>% 
        add_linear_penalties(penalty = 5, data = my.tracks) %>% #prefer cells with real trajectory data in current season
        add_linear_penalties(penalty = 2.5, data = any.tracks) %>% #prefer cells cells with real trajectory data in any season
        #add_linear_penalties(penalty = 5, data = binary.ves) %>% #prefer cells which are not regularly used routes across any season
        add_boundary_penalties(penalty = my.boundary.penalty, edge_factor = 0) %>% 
        add_gurobi_solver(time_limit = 1000, threads = 8, gap = 0.1) %>% 
        solve()
      
      names(priorit.solve) <- paste0("priorit.", i)
      names(priorit.solve.free) <- paste0("priorit.", i)
      
      assign(paste0("priorit.solve.", covi), priorit.solve, environment())
      assign(paste0("priorit.solve.free.", covi), priorit.solve.free, environment())
      
    }
    
    
    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####
    
    my.params <- raster(my.tracks) # store these in first entry in empty raster
    my.params[my.params][1] <- round(my.boundary.penalty, 4)
    
    names(my.ves.ras) <- "ves.dens"
    names(my.spat.pred) <- "spat.pred"
    names(my.tracks) <- "my.tracks"
    names(binary.ves) <- "binary.ves"
    names(any.tracks) <- "any.tracks"
    names(my.params) <- "my.boundary.penalty"
    
    my.priority.zones.1[[length(my.priority.zones.1)+1]] <- stack(priorit.solve.1, my.ves.ras, my.spat.pred, my.tracks, binary.ves, any.tracks, my.params)
    my.priority.zones.2[[length(my.priority.zones.2)+1]] <- stack(priorit.solve.2, my.ves.ras, my.spat.pred, my.tracks, binary.ves, any.tracks, my.params)
    my.priority.zones.free.1[[length(my.priority.zones.free.1)+1]] <- stack(priorit.solve.free.1, my.ves.ras, my.spat.pred, my.tracks, binary.ves, any.tracks, my.params)
    my.priority.zones.free.2[[length(my.priority.zones.free.2)+1]] <- stack(priorit.solve.free.2, my.ves.ras, my.spat.pred, my.tracks, binary.ves, any.tracks, my.params)
    
  }
  
  all.priority.zones <- list(my.priority.zones.1, my.priority.zones.2, my.priority.zones.free.1, my.priority.zones.free.2)
  names(all.priority.zones) <- c("cost.zones.1", "cost.zones.2", "free.zones.1", "free.zones.2")
  
  return(all.priority.zones)
  
}

plot(spat.pred$iseason.ras, spat.pred$idate) # iseason ras is now same range as index.idate/idate

my.zones <- create.priority.zones(spatial.pred = spat.pred, vessel.data = vessel.data.all, 
                                  max.clustering = 0.5, cover = c(0.10, 0.30), # onyl works for 2 covers atm 
                                  resolution.fact = resolution.fact)

my.zones$cost.zones.1[[1]]
my.zones$free.zones.1[[1]]


how.much.protect <- function(zones){
  
  how.much.all <- list()
  
  for( i in 1:length(zones)){
    
    my.zone <- zones[[i]]
    
    my.zone <- lapply(my.zone, function(x){
      mpa.summary.stats <- function(data){
        
        protect.index <- which(data[[1]][] == 1)
        
        how.much.protc <- data.frame(which.zone = names(zones)[i],
                                     iseason.ras = substr(names(data)[1], 9, 10),
                                     how.much.area = length(which((data[[1]][] == 1))) / length(which(!is.na(data[[1]][]))),
                                     how.much.spat = (sum(data[["spat.pred"]][protect.index], na.rm = T) / sum(data[["spat.pred"]][], na.rm = T)),
                                     #how.much.idv = (sum((data[["my.tracks"]]-1)[protect.index], na.rm = T) / sum((data[["my.tracks"]]-1)[], na.rm = T)),
                                     #conv mytracks because its inverse in calc
                                     #removed idv  because i dont think its itneresting
                                     how.much.ves = (sum(data[["ves.dens"]][protect.index], na.rm = T) / sum(data[["ves.dens"]][], na.rm = T)))
        
        return(how.much.protc)
      }
      return(mpa.summary.stats(x))
    })
    my.zone <- do.call("rbind", my.zone)
    how.much.all[[(length(how.much.all) + 1)]] <- my.zone
  }
  
  return(how.much.all)
  
}

how.much <- do.call("rbind", how.much.protect(my.zones))

how.much.summary <- data.frame(which.zone = c("cost.zones.1", "cost.zones.1", "free.zones.1", "free.zones.2"),
                               mean.area.cover = c(mean(subset(how.much, which.zone == "cost.zones.1")$how.much.area),
                                                   mean(subset(how.much, which.zone == "cost.zones.2")$how.much.area),
                                                   mean(subset(how.much, which.zone == "free.zones.1")$how.much.area),
                                                   mean(subset(how.much, which.zone == "free.zones.2")$how.much.area)),
                               mean.spat.cover = c(mean(subset(how.much, which.zone == "cost.zones.1")$how.much.spat),
                                                   mean(subset(how.much, which.zone == "cost.zones.2")$how.much.spat),
                                                   mean(subset(how.much, which.zone == "free.zones.1")$how.much.spat),
                                                   mean(subset(how.much, which.zone == "free.zones.2")$how.much.spat)),
                               mean.ves.cover = c(mean(subset(how.much, which.zone == "cost.zones.1")$how.much.ves),
                                                  mean(subset(how.much, which.zone == "cost.zones.2")$how.much.ves),
                                                  mean(subset(how.much, which.zone == "free.zones.1")$how.much.ves),
                                                  mean(subset(how.much, which.zone == "free.zones.2")$how.much.ves)),
                               
                               sd.area.cover = c(sd(subset(how.much, which.zone == "cost.zones.1")$how.much.area),
                                                 sd(subset(how.much, which.zone == "cost.zones.2")$how.much.area),
                                                 sd(subset(how.much, which.zone == "free.zones.1")$how.much.area),
                                                 sd(subset(how.much, which.zone == "free.zones.2")$how.much.area)),
                               sd.spat.cover = c(sd(subset(how.much, which.zone == "cost.zones.1")$how.much.spat),
                                                 sd(subset(how.much, which.zone == "cost.zones.2")$how.much.spat),
                                                 sd(subset(how.much, which.zone == "free.zones.1")$how.much.spat),
                                                 sd(subset(how.much, which.zone == "free.zones.2")$how.much.spat)),
                               sd.ves.cover = c(sd(subset(how.much, which.zone == "cost.zones.1")$how.much.ves),
                                                sd(subset(how.much, which.zone == "cost.zones.2")$how.much.ves),
                                                sd(subset(how.much, which.zone == "free.zones.1")$how.much.ves),
                                                sd(subset(how.much, which.zone == "free.zones.2")$how.much.ves)))

how.much <- #make long format
  split(how.much, how.much$which.zone) %>%
  lapply(., function(x){
    x <- data.frame(which.zone = rep(x$which.zone, 3),
                    iseason.ras = rep(x$iseason.ras, 3),
                    what = rep(c("area", "spat", "ves"), each = 12),
                    how.much = c(x$how.much.area, x$how.much.spat, x$how.much.ves))
    
    return(x)
  }) %>%
  do.call("rbind", .)
row.names(how.much) <- NULL

order <- 1:12
for( i in unique(how.much$which.zone)){
  my.plot <- ggplot(data = subset(how.much, which.zone == i)) +
    geom_line(lwd = 2, stat = "identity",
              aes(x = iseason.ras, y = how.much, group = what, colour = what)) +
    scale_x_discrete(limits = function(x){order[order %in% x]}) +
    scale_y_continuous(limits = c(0, max(how.much$how.much)))
  assign(paste0("plot.", i), my.plot, envir = globalenv())
  rm(my.plot)
}


dates.match <- 
  unique(fit.mpm[,c("iseason.ras", "index.idate", "date")]@data) %>%
  split(., .$iseason.ras) %>%
  lapply(., function(x){
    return(data.frame(iseason.ras = unique(x$iseason.ras),
                      index.idate = unique(x$index.idate),
                      date = median(as.Date(paste0("1966-", substr(x$date,6,10)), tryFormats = "%Y-%m-%d"))
    ))
  }) %>%
  do.call("rbind", .)

#notesinfunction
summary.output.paplan <- lapply(season.split, function(x){
  
  #x <- season.split[[6]]
  
  iseas <- unique(x$ISeason)
  print(iseas)
  
  how.many.km <- sum(x$dis.trav)
  how.many.hr <- sum(x$time.trav)
  
  cost.zones.1 <- my.zones$cost.zones.1[[iseas]][[paste0("priorit.", iseas)]]
  cost.zones.2 <- my.zones$cost.zones.2[[iseas]][[paste0("priorit.", iseas)]]
  
  free.zones.1 <- my.zones$free.zones.1[[iseas]][[paste0("priorit.", iseas)]]
  free.zones.2 <- my.zones$free.zones.2[[iseas]][[paste0("priorit.", iseas)]]
  
  names(cost.zones.1) <- "priorit"
  names(cost.zones.2) <- "priorit"
  names(free.zones.1) <- "priorit"
  names(free.zones.2) <- "priorit"
  
  inside.the.lines <- function(data){
    inside <- 
      rasterToPolygons(data, dissolve = TRUE) %>%
      subset(., priorit == 1) %>%
      st_as_sf(.) %>%
      st_intersection(st_as_sf(x), .)
    return(inside)
    
    #outside.cost.1 <- #not that useful because of data values being for full line segement, ki.e. there not recalcualted
    #  rasterToPolygons(cost.zones.1, dissolve = TRUE) %>%
    #  subset(., priorit == 0) %>%
    #  st_as_sf(.) %>%
    #  st_intersection(st_as_sf(x), .)
  }
  
  cost.zones.1 <- inside.the.lines(cost.zones.1)
  cost.zones.2 <- inside.the.lines(cost.zones.2)
  free.zones.1 <- inside.the.lines(free.zones.1)
  free.zones.2 <- inside.the.lines(free.zones.2)
  
  
  
  total.dis <- sum(st_length(st_as_sf(x))) # total dis travelled
  total.dis.in.cost.1 <- sum(st_length(cost.zones.1)) # total dis travelled inside protected area
  total.dis.in.cost.2 <- sum(st_length(cost.zones.2)) # total dis travelled inside protected area
  total.dis.in.free.1 <- sum(st_length(free.zones.1)) # total dis travelled inside protected area
  total.dis.in.free.2 <- sum(st_length(free.zones.2)) # total dis travelled inside protected area
  
  
  summary.output <- data.frame(iseason = rep(unique(x$ISeason), times = 5), 
                               mid.date = dates.match$date[which(dates.match$iseason.ras == unique(x$ISeason))],
                               zone = c("total", 
                                        "cost.1", 
                                        "cost.2", 
                                        "free.1", 
                                        "free.2"),
                               total.dis = as.numeric(c(total.dis,
                                                        total.dis.in.cost.1,
                                                        total.dis.in.cost.2,
                                                        total.dis.in.free.1,
                                                        total.dis.in.free.2)),
                               perc.dis = as.numeric(c(1,
                                                       total.dis.in.cost.1 / total.dis,
                                                       total.dis.in.cost.2 / total.dis,
                                                       total.dis.in.free.1 / total.dis,
                                                       total.dis.in.free.2 / total.dis)),
                               mean.speed.kmh = as.numeric(c(mean(x$Calc.Speed.kmh),
                                                             mean(cost.zones.1$Calc.Speed.kmh),
                                                             mean(cost.zones.2$Calc.Speed.kmh),
                                                             mean(free.zones.1$Calc.Speed.kmh), 
                                                             mean(free.zones.2$Calc.Speed.kmh))),
                               sd.speed = as.numeric(c(sd(x$Calc.Speed.kmh),
                                                       sd(cost.zones.1$Calc.Speed.kmh),
                                                       sd(cost.zones.2$Calc.Speed.kmh),
                                                       sd(free.zones.1$Calc.Speed.kmh), 
                                                       sd(free.zones.2$Calc.Speed.kmh))))
  
  return(summary.output)
})
summary.output.paplan <- do.call("rbind", summary.output.paplan)
summary.output.paplan$cost.free <- substr(summary.output.paplan$zone, 1, 4)
summary.output.paplan$sum.cover <- substr(summary.output.paplan$zone, 6, 7)
summary.output.paplan$cost.free[which(summary.output.paplan$zone == "total")] <- "total"
summary.output.paplan$sum.cover[which(summary.output.paplan$zone == "total")] <- 3

##################################################
##need continuous dates for nice looking plot
summary.output.paplan.forplot <- 
  data.frame(date = seq(from = as.Date(paste0("1966-", 
                                              substr(unique(fit.mpm$date[which(fit.mpm$idate == 1)]), 6, 10)), 
                                       tryFormats = "%Y-%m-%d"),
                        to = as.Date(paste0("1967-", 
                                            substr(unique(fit.mpm$date[which(fit.mpm$idate == 1)]), 6, 10)), 
                                     tryFormats = "%Y-%m-%d")-1, by = 1),
             idate = 1:365)
summary.output.paplan.forplot$iseason.ras <- inla.group(summary.output.paplan.forplot$idate, n = 12, idx.only = TRUE) 
# same as iseason# see hpc


summary.output.paplan.forplot.ls <- list()
for( i in 1:length(unique(summary.output.paplan$zone))){
  summary.output.paplan.forplot.ls[[i]] <- summary.output.paplan.forplot
  summary.output.paplan.forplot.ls[[i]]$zone <- unique(summary.output.paplan$zone)[i]
}
summary.output.paplan.forplot <- summary.output.paplan.forplot.ls
rm(summary.output.paplan.forplot.ls)

summary.output.paplan.forplot <- 
  lapply(summary.output.paplan.forplot, function(x){
    
    my.targ <- subset(summary.output.paplan, zone == unique(x$zone))
    my.targ$zone <- NULL
    
    x <- merge(x, my.targ, by.x = "iseason.ras", by.y = "iseason")
    return(x)
  }) %>%
  do.call("rbind", .)

#make all same year so can plot on x scale
summary.output.paplan.forplot$plot.date <- as.Date(paste0("1966", substr(summary.output.paplan.forplot$date, 5, 10)),
                                                   tryFormats = "%Y-%m-%d")

#save.image("Telemetry.spatial.predictions.protected.areas.RData")

