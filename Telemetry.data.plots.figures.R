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

load("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Telemetry.spatial.predictions.protected.areas.RData")

################################################################################################################################################
#this script 
#- creates plots and figures

#################################################################################################

river.col <- "deepskyblue3"
continental.col <- "tan2"
haul.col <- "mediumseagreen"
sst.col <- "coral3"

pals <- data.frame(river.col = river.col,
                   continental.col = continental.col,
                   haul.col = haul.col,
                   sst.col = sst.col)

#posterior marginal for fixed effects plots

theme.standards <- 
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11 ),#, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12),
        plot.title = element_text(size=12, vjust = -0.15),
        #axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,0.9,0), "cm")) 

########################################################################################################################################

model = model.time
covs = names(model.time$marginals.fixed)
names = c("Distance from\nriver mouths",
          "Distance from\ncontinental shelf",
          #"Distance from haulout",
          "Sea surface\ntemperature")
my.pal = c(pals$river.col, pals$continental.col, pals$sst.col)
rm(model, covs, names, my.pal)


create.fixed.plot <- function(model, covs, names, my.pal){
  ###################################################################
  fixed.pred <- 
    model$marginals.fixed %>% 
    lapply(., function(y){
      return(data.frame(inla.smarginal(y)))
    })
  for(i in 1:length(fixed.pred)){
    fixed.pred[[i]]$effect <-names(fixed.pred)[i]
  }
  fixed.pred <- do.call("rbind", fixed.pred)
  ###################################################################
  my.quantiles <- 
    model$marginals.fixed %>% 
    lapply(., function(y){
      return(data.frame(x0.025 = as.numeric(inla.qmarginal(0.025, y)), 
                        x0.25 = as.numeric(inla.qmarginal(0.25, y)),
                        x0.50 = as.numeric(inla.qmarginal(0.5, y)),
                        x0.75 = as.numeric(inla.qmarginal(0.75, y)),
                        x0.975 = as.numeric(inla.qmarginal(0.975, y))))
    }) 
  for(i in 1:length(my.quantiles)){
    my.quantiles[[i]]$effect <- names(my.quantiles)[i]
  }
  my.quantiles <- do.call("rbind", my.quantiles)
  ###################################################################
  
  
  fixed.plot <- ggplot() +
    geom_hline(alpha = 1, yintercept = 0 , linetype = "dashed", lwd = 1.2, colour = "grey68") +
    see::geom_violinhalf(data = fixed.pred, width = 2, colour = "gray98", scale = "area", 
                         alpha = 0.75, adjust = 1, trim=TRUE, aes(x=effect, y=x, group = effect, fill = effect)) +
    scale_y_continuous(limits = c(-max(abs(fixed.pred$x) * 1.25), max(abs(fixed.pred$x) * 1.25))) +
    scale_x_discrete(limits = covs,
                     labels = names) +
    scale_fill_manual(limits = covs,
                      values = my.pal) +
    geom_errorbar(data = my.quantiles, width = 0.08, lwd = 1.2, colour = "gray40",
                  aes(x = effect, group = effect, ymin = x0.025, ymax = x0.975)) + 
    geom_errorbar(data = my.quantiles, width = NA, lwd = 5, colour = "gray80",
                  aes(x = effect, group = effect, ymin = x0.25, ymax = x0.75)) +
    #geom_errorbar(data = my.quantiles, width = NA, lwd = 4, colour = "white",
    #              aes(x = effect, group = effect, ymin = x0.50-1e-09, ymax = x0.50+1e-09)) +
    #geom_point(data = my.quantiles, colour = "gray99", size = 2,
    #           aes(x = effect, group = effect, y = x0.50)) +
    ylab(label = "Coefficient") +
    coord_flip() +
    labs(title = "Posterior marginal distribution for the\nfixed effects of environmental covariates") +
    theme.standards +
    theme(axis.text.y= element_text(hjust=0.5, size = 11, angle = 90, vjust = 0),
          axis.title.y=element_blank(),
          legend.position = "none")
  return(fixed.plot)
}

fixed.plot <- create.fixed.plot(model = model.time,
                                covs = names(model.time$marginals.fixed),
                                names = c("Distance from\nriver mouths",
                                          "Distance from\ncontinental shelf",
                                          #"Distance from haulout",
                                          "Sea surface\ntemperature"),
                                my.pal = c(pals$river.col, pals$continental.col, pals$sst.col))

fixed.plot


########################################################################################################################################

create.temporal.plot <- function(model, full.or.partial, my.pal) {
  
  date.pred <- data.frame(idate = seq(from = min(fit.mpm$idate), to = max(fit.mpm$idate), by = 1),
                          date = seq(from = as.Date(paste0("1966-", unique(substr(fit.mpm$date[which(fit.mpm$idate == 1)], 6, 10))), 
                                                    tryFormats = "%Y-%m-%d"),
                                     to = as.Date(paste0("1967-", unique(substr(fit.mpm$date[which(fit.mpm$idate == max(fit.mpm$idate))], 6, 10))), 
                                                  tryFormats = "%Y-%m-%d"),
                                     by = 1,),
                          cent.sq.riv = 1,
                          cent.sq.dis.con = 1,
                          cent.sq.dis.haul.any = 1, # 1 so can estimate full plus diff,
                          cent.sq.dis.haul.any = 1, # 1 so can estimate full plus diff
                          cent.sst = 1, 
                          bath.abs = 1
  )
  
  date.pred$plot.date <- paste0("1996-", substr(date.pred$date, 6, 10))
  date.pred <- date.pred[order(date.pred$plot.date),]
  date.pred$plot.idate <- 1:nrow(date.pred)
  row.names(date.pred) <- NULL
  date.pred$date <- substr(date.pred$date, 6, 10)
  date.pred$plot.date <- NULL
  
  #date.pred.smooth <- predict(model, date.pred, ~ season.smooth, include = "season.smooth")
  
  if(full.or.partial == "full"){
    
    date.pred.riv <- predict(model, date.pred, ~ riv.cov + riv.cov.time, include = c("riv.cov", "riv.cov.time"))
    date.pred.con <- predict(model, date.pred, ~ con.cov + con.cov.time, include = c("con.cov", "con.cov.time"))
    #date.pred.haul <- predict(model, date.pred, ~ haul.cov + haul.cov.time, include = c("haul.cov", "haul.cov.time"))
    
  }
  
  if(full.or.partial == "partial"){
    
    date.pred.riv <- predict(model, date.pred, ~ riv.cov.time, include = c("riv.cov.time"))
    date.pred.con <- predict(model, date.pred, ~ con.cov.time, include = c("con.cov.time"))
    #date.pred.haul <- predict(model, date.pred, ~ haul.cov.time, include = c("haul.cov.time"))
    
  }
  
  pred.names <- c("date.pred.riv", "date.pred.con") #, "date.pred.haul") # "date.pred.smooth", 
  
  date.pred <- list()
  for( i in 1:length(pred.names)){
    date.pred[[i]] <-
      ggplot() +
      geom_ribbon(data = eval(parse(text = pred.names[i])), alpha = 0.75, fill = my.pal[i],
                  aes(x = plot.idate, ymin = q0.025, ymax = q0.975)) +
      geom_line(data = eval(parse(text = pred.names[i])), colour = "grey40", aes(x = plot.idate, y = mean)) + 
      scale_x_continuous(breaks=which(substr(date.pred.riv$date, 4, 5) == "01"),
                         labels = c("January", "February", "March", "April", "May", "June", "July", "August", 
                                    "September", "October", "November", "December")) +
      ylab(label = "Coefficient")
    
    if(pred.names[i] != "date.pred.smooth"){
      date.pred[[i]] <- date.pred[[i]] + 
        geom_hline(yintercept = 0 ,linetype = "dashed", lwd = 1.2, colour = "grey40") 
    }
    if(pred.names[i] == "date.pred.smooth"){
      date.pred[[i]] <- date.pred[[i]] + 
        labs(title = "Seasonal trend") +
        theme.standards +
        theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust=1),
              axis.title.x=element_blank(),
              plot.margin = unit(c(0,0,0.5,0.5), "cm"))
    }
    if(grepl("riv", pred.names[i])){
      date.pred[[i]] <- date.pred[[i]] + 
        labs(title = "Distance from river") +
        theme.standards +
        theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust=1),
              axis.title.x=element_blank(),
              plot.margin = unit(c(0,0,0.5,0.5), "cm"))
    }
    if(grepl("con", pred.names[i])){
      date.pred[[i]] <- date.pred[[i]] + 
        labs(title = "Distance from continental shelf") +
        theme.standards +
        theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust=1),
              axis.title.x=element_blank(),
              plot.margin = unit(c(0,0,0.5,0.5), "cm"))
      
    }
    if(grepl("haul", pred.names[i])){
      date.pred[[i]] <- date.pred[[i]] + 
        labs(title = "Distance from haulout") +
        theme.standards +
        theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust=1),
              axis.title.x=element_blank(),
              plot.margin = unit(c(0,0,0.5,0.5), "cm"))
    }
  }
  
  temporal.plot <- cowplot::plot_grid(plotlist = date.pred[length(date.pred):1], ncol = 1)
  
  return(temporal.plot)
}

temporal.plot.full <- create.temporal.plot(model.time, 
                                           full.or.partial = "full",
                                           my.pal = c(pals$river.col, pals$continental.col))


mpa.dis.perc.cover.plots <- 
  ggplot(data = subset(summary.output.paplan.forplot, zone != "total")) +
  geom_line(lwd = 2, 
            aes(x = plot.date, y = perc.dis*100, group = zone, colour = sum.cover, linetype = cost.free)) +
  scale_x_continuous(breaks = seq(from = as.Date("1966-01-01", tryFormats = "%Y-%m-%d"),
                                  to = as.Date( "1966-12-31", tryFormats = "%Y-%m-%d"), length.out = 12),
                     labels = c("January", "February", "March", "April", "May", "June", "July", "August", 
                                "September", "October", "November", "December")) +
  scale_colour_discrete(name = "MPA plan",
                        labels = c("10%", "30%"), 
                        type = c("#F0A500", "#05595B")) + 
  scale_linetype_discrete(name = element_blank(),
                          labels = c("Vessel\npenalized", "Unpenalized")) + 
  theme.standards +
  labs(title = "Percentage of total distance travelled\nwithin each MPA plan") +
  ylab(label = "Percentage of total distance travelled") +
  theme(axis.text.x= element_text(size = 11, angle = 45, vjust = 0.5),
        axis.title.x=element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(order = 1), 
         linetype = guide_legend(order = 2))

mpa.speed.plots <- 
  ggplot(data = subset(summary.output.paplan.forplot, zone != "total")) +
  geom_line(lwd = 2, 
            aes(x = plot.date, y = mean.speed.kmh, group = zone, colour = sum.cover, linetype = cost.free)) +
  scale_x_continuous(breaks = seq(from = as.Date("1966-01-01", tryFormats = "%Y-%m-%d"),
                                  to = as.Date( "1966-12-31", tryFormats = "%Y-%m-%d"), length.out = 12),
                     labels = c("January", "February", "March", "April", "May", "June", "July", "August", 
                                "September", "October", "November", "December")) +
  scale_colour_discrete(name = "MPA plan",
                        labels = c("10%", "30%"), 
                        type = c("#F0A500", "#05595B")) + 
  scale_linetype_discrete(name = element_blank(),
                          labels = c("Vessel\npenalized", "Unpenalized")) + 
  theme.standards +
  labs(title = "Mean vessel speed within each MPA plan\n ") +
  ylab(label = "Mean vessel speed (kmh)") +
  theme(axis.text.x= element_text(size = 11, angle = 45, vjust = 0.5),
        axis.title.x=element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) 


###################

#create spatial maps
spat.pal <- inlmisc::GetColors(200, scheme = "smooth rainbow", reverse = FALSE, stops = c(0.2,0.85))

theme.spatials <- theme(legend.position = "none",
                        plot.margin = unit(c(0,0,0,0), "cm"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        panel.grid=element_blank(), 
                        panel.background=element_rect(fill = "transparent",colour = NA),
                        panel.border=element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        plot.title = element_text(size = 11))


model.maps <- list()
heat.maps <- list()
vessel.maps <- list()
priorit.maps.1 <- list()
priorit.maps.2 <- list()
for( y in 1:max(fit.mpm$iseason.ras)){
  
  my.ves.pix <- as(mask(sqrt(vessel.data.all[[y]]), caspian.shape), "SpatialPixelsDataFrame")
  
  dates <- unique(fit.mpm[,c("date", "iseason.ras")]@data)
  dates <- subset(dates, iseason.ras == y)
  dates$date <- as.Date(paste0("1966", substr(dates$date, 5, 10)))
  
  if ((12 %in% as.numeric(substr(dates$date, 6, 7))) & (1 %in% as.numeric(substr(dates$date, 6, 7)))){
    title <- paste0(substr(min(dates$date[which(substr(dates$date, 6, 7) == "12")]), 6, 10), " to\n", 
                    substr(max(dates$date[which(substr(dates$date, 6, 7) == "01")]), 6, 10))
  } else {
    title <- paste0(min(substr(dates$date, 6, 10)), " to\n", max(substr(dates$date, 6, 10)))
  } # if loops over dec and jan, need to change title order
  
  my.zone.1 <- my.zones$cost.zones.1[[y]][[1]]
  my.zone.free.1 <- my.zones$free.zones.1[[y]][[1]]
  my.zone.2 <- my.zones$cost.zones.2[[y]][[1]]
  my.zone.free.2 <- my.zones$free.zones.2[[y]][[1]]
  
  smoothen.mpas <- function(x, res){
    names(x) <- "priorit"
    x <- rasterToPolygons(x, dissolve = TRUE) %>%
      subset(., priorit == 1) %>% 
      smoothr::smooth(., method = "ksmooth", smoothness = 5) %>%
      gBuffer(x, byid = T, width = res)
    #smoothen out MPA
    return(x)
  }
  my.zone.1 <- smoothen.mpas(my.zone.1, res = res(my.zones$cost.zones.1[[y]][[1]])[1]/2)
  my.zone.free.1 <- smoothen.mpas(my.zone.free.1, res = res(my.zones$cost.zones.1[[y]][[1]])[1]/2)
  my.zone.2 <- smoothen.mpas(my.zone.2, res = res(my.zones$cost.zones.1[[y]][[1]])[1]/2)
  my.zone.free.2 <- smoothen.mpas(my.zone.free.2, res = res(my.zones$cost.zones.1[[y]][[1]])[1]/2)
  
  my.zone.1$free <- 0
  my.zone.free.1$free <- 1
  my.zone.2$free <- 0
  my.zone.free.2$free <- 1
  
  my.zone.1 <- rbind(my.zone.1, my.zone.free.1)
  my.zone.2 <- rbind(my.zone.2, my.zone.free.2)
  
  rm(my.zone.free.1, my.zone.free.2)
  
  my.zone.1 <- disaggregate(my.zone.1)
  my.zone.1$id <- 1:nrow(my.zone.1)
  my.zone.1 <- crop(my.zone.1, caspian.shape)
  
  my.zone.2 <- disaggregate(my.zone.2)
  my.zone.2$id <- 1:nrow(my.zone.2)
  my.zone.2 <- crop(my.zone.2, caspian.shape)
  
  #################################################
  
  heat <- aggregate(subset(fit.mpm, iseason.ras == y)[,"g.auc.stan.inv"], 
                    subset(spat.pred[,"iseason.ras"], iseason.ras == y), 
                    FUN = function(x){
                      return(mean(x, na.rm = T))})
  
  
  spat.pred@grid@cellsize
  heat@grid@cellsize
  my.ves.pix@grid@cellsize
  
  
  model.map <- 
    ggplot() +
    gg(subset(spat.pred, iseason.ras == y), aes(fill = median)) +
    scale_fill_gradientn(limits = c(0,1), colours = spat.pal,
                         name = "Predicted\n1 - g") +
    #labs(title = paste0("Model map for ", paste(range(sort(substr(unique(subset(fit.mpm, iseason.ras == y)$date), 6, 10))), 
    #                    collapse = " to "))) +
    #scale_x_continuous(range(spat.pred@coords[,1])) +
    #scale_y_continuous(range(spat.pred@coords[,2])) +
    coord_cartesian(xlim = range(spat.pred@coords[,1]), ylim = range(spat.pred@coords[,2])) +
    labs(title = title) +
    theme.standards +
    theme.spatials 
  
  
  heat.map <- 
    ggplot() + 
    gg(heat, aes(fill = g.auc.stan.inv)) +
    scale_fill_gradientn(limits = c(0,1), colours = spat.pal,
                         name = "Observed\n1 - g") +
    #labs(title = paste0("Heat map for ", paste(range(sort(substr(unique(subset(fit.mpm, iseason.ras == y)$date), 6, 10))), 
    #                                           collapse = " to "))) +
    #scale_x_continuous(range(spat.pred@coords[,1])) +
    #scale_y_continuous(range(spat.pred@coords[,2])) +
    coord_cartesian(xlim = range(spat.pred@coords[,1]), ylim = range(spat.pred@coords[,2])) +
    labs(title = "") +
    theme.standards +
    theme.spatials
  
  vessel.map <- 
    ggplot() +
    gg(my.ves.pix) +
    #gg(as(mask(vessel.data.all[[y]], caspian.shape), "SpatialPixelsDataFrame")) +
    #gg(as(mask(vessel.data.all[[y]], caspian.shape), "SpatialPixelsDataFrame")) +
    scale_fill_gradientn(colours = spat.pal,
                         limits = c(0,max(sqrt(vessel.data.all[]))),
                         name = "Sqrt vessels per km2") +
    #labs(title = paste0("Model map for ", paste(range(sort(substr(unique(subset(fit.mpm, iseason.ras == y)$date), 6, 10))), 
    #                    collapse = " to "))) +
    coord_cartesian(xlim = range(spat.pred@coords[,1]), ylim = range(spat.pred@coords[,2])) +
    labs(title = "") +
    theme.standards +
    theme.spatials
  
  
  priorit.map.1 <- 
    ggplot() +
    gg(caspian.shape, fill = NA, colour = NA) +
    gg(heat, fill = "grey50") +
    gg(data = my.zone.1, alpha= 0.4, lwd = 0.25, aes(fill = as.factor(free), colour = as.factor(free))) + #, x = x, y = y)) +
    scale_fill_manual(limits = c("0","1"),
                      labels = c("Vessel\npenalized", "Unpenalized"),
                      values = c("#FF4C29", "#035397"),
                      name = "10% MPA plan") +
    scale_colour_manual(limits = c("0","1"),
                        values = c("#FF4C29", "#035397")) +
    coord_cartesian(xlim = range(spat.pred@coords[,1]), ylim = range(spat.pred@coords[,2])) +
    labs(title = "") +
    theme.standards +
    theme.spatials
  
  priorit.map.2 <- 
    ggplot() +
    gg(caspian.shape, fill = NA, colour = NA) +
    gg(heat, fill = "grey50") +
    gg(data = my.zone.2, alpha= 0.4, lwd = 0.25, aes(fill = as.factor(free), colour = as.factor(free))) + #, x = x, y = y)) +
    scale_fill_manual(limits = c("0","1"),
                      labels = c("Vessel\npenalized", "Unpenalized"),
                      values = c("#FF4C29", "#035397"),
                      name = "30% MPA plan") +
    scale_colour_manual(limits = c("0","1"),
                        values = c("#FF4C29", "#035397")) +
    coord_cartesian(xlim = range(spat.pred@coords[,1]), ylim = range(spat.pred@coords[,2])) +
    labs(title = "") +
    theme.standards +
    theme.spatials
  
  model.leg <- cowplot::get_legend(model.map + theme(legend.position = "right"))
  heat.leg <- cowplot::get_legend(heat.map + theme(legend.position = "right"))
  vessel.leg <- cowplot::get_legend(vessel.map + theme(legend.position = "right"))
  priorit.leg.1 <- cowplot::get_legend(priorit.map.1 + theme(legend.position = "right"))
  priorit.leg.2 <- cowplot::get_legend(priorit.map.2 + theme(legend.position = "right"))
  
  model.maps[[y]] <- model.map
  heat.maps[[y]] <- heat.map
  vessel.maps[[y]] <- vessel.map
  priorit.maps.1[[y]] <- priorit.map.1
  priorit.maps.2[[y]] <- priorit.map.2
  
  #spatial.maps[[y]] <- cowplot::plot_grid(model.map, heat.map, vessel.map, priorit.map, nrow = 1)
  
} # 1:max(fit.mpm$iseason.ras)

model.maps <- model.maps[c(9:12,1:8)] # reorder to make more sense of date range, i.e. starting at january ish
heat.maps <- heat.maps[c(9:12,1:8)]
priorit.maps.1 <- priorit.maps.1[c(9:12,1:8)]
priorit.maps.2 <- priorit.maps.2[c(9:12,1:8)]

model.maps[[length(model.maps)+1]] <- model.leg
heat.maps[[length(heat.maps)+1]] <- heat.leg
priorit.maps.1[[length(priorit.maps.1)+1]] <- priorit.leg.1
priorit.maps.2[[length(priorit.maps.2)+1]] <- priorit.leg.2




############################################################################################################
#make season averaged SST 

season.dates <- unique(data.frame(index.idate = fit.mpm$index.idate, date = fit.mpm$date))
season.dates$date <- gsub("-", "", season.dates$date)

SST.season.mean <- list()
for( i in sort(unique(fit.mpm$index.idate))){
  SST.season.mean[[i]] <- 
    SST[which(substr(SST, 56, 63) %in% unique(subset(season.dates, index.idate == i)$date))] %>%
    stack(.) %>%
    mean(., na.rm = T) %>%
    projectRaster(from = ., to = distance.any.ras) - 273.15 # celcius
}
SST.season.mean <- 
  stack(SST.season.mean) %>%
  as(., "SpatialPixelsDataFrame") %>%
  crop(., caspian.shape)
names(SST.season.mean) <- gsub("layer", "mean.temp", names(SST.season.mean))

SST.season.sd <- list()
for( i in sort(unique(fit.mpm$index.idate))){
  SST.season.sd[[i]] <- 
    SST[which(substr(SST, 56, 63) %in% unique(subset(season.dates, index.idate == i)$date))] %>%
    stack(.) %>%
    calc(., fun = sd, na.rm = T) %>%
    projectRaster(from = ., to = distance.any.ras) 
}
SST.season.sd <- 
  stack(SST.season.sd) %>%
  as(., "SpatialPixelsDataFrame") %>%
  crop(., caspian.shape)
names(SST.season.sd) <- gsub("layer", "sd.temp", names(SST.season.sd))

SST.plots <- as.list(sort(unique(fit.mpm$index.idate)))
SST.plots <- lapply(SST.plots, function(x){
  for.plot.m <- SST.season.mean
  for.plot.sd <- SST.season.sd
  
  names(for.plot.m)[which(names(for.plot.m) == paste0("mean.temp.", x))] <- "mean.for.plot"
  names(for.plot.sd)[which(names(for.plot.sd) == paste0("sd.temp.", x))] <- "sd.for.plot"
  
  x.m <- 
    ggplot() +
    gg(for.plot.m, aes(fill = mean.for.plot)) +
    scale_fill_gradientn(colours = spat.pal, name = "Sea surface\ntemperature (°C)",
                         limits = c(0, max(SST.season.mean$mean.temp.1, SST.season.mean$mean.temp.2, SST.season.mean$mean.temp.3, SST.season.mean$mean.temp.4, SST.season.mean$mean.temp.5,
                                           SST.season.mean$mean.temp.6, SST.season.mean$mean.temp.7, SST.season.mean$mean.temp.8, SST.season.mean$mean.temp.9, SST.season.mean$mean.temp.10,
                                           SST.season.mean$mean.temp.11))) +
    theme.standards +
    theme(legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  x.sd <- 
    ggplot() +
    gg(for.plot.sd, aes(fill = sd.for.plot)) +
    scale_fill_gradientn(colours = spat.pal, name = "Sea surface\ntemperature (°C)",
                         limits = c(min(SST.season.sd$sd.temp.1, SST.season.sd$sd.temp.2, SST.season.sd$sd.temp.3, SST.season.sd$sd.temp.4, SST.season.sd$sd.temp.5,
                                        SST.season.sd$sd.temp.6, SST.season.sd$sd.temp.7, SST.season.sd$sd.temp.8, SST.season.sd$sd.temp.9, SST.season.sd$sd.temp.10,
                                        SST.season.sd$sd.temp.11), 
                                    max(SST.season.sd$sd.temp.1, SST.season.sd$sd.temp.2, SST.season.sd$sd.temp.3, SST.season.sd$sd.temp.4, SST.season.sd$sd.temp.5,
                                        SST.season.sd$sd.temp.6, SST.season.sd$sd.temp.7, SST.season.sd$sd.temp.8, SST.season.sd$sd.temp.9, SST.season.sd$sd.temp.10,
                                        SST.season.sd$sd.temp.11))) +
    theme.standards +
    theme(legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  
  
  return(cowplot::plot_grid(x.m, x.sd, ncol = 1))
})

cowplot::plot_grid(plotlist = SST.plots, nrow = 1)

date.titles <- as.list(1:max(fit.mpm$index.idate))
date.titles <- lapply(date.titles, function(x){
  dates <- subset(season.dates, index.idate == x)$date
  
  min.d <- as.Date(paste0(1966, substr(dates, 5, 10)), tryFormats = c("%Y%m%d"))
  max.d <- as.Date(paste0(1967, substr(dates, 5, 10)), tryFormats = c("%Y%m%d"))
  
  if((sum(grepl("12", substr(max.d, 6, 7)) != 0)) & (sum(grepl("01", substr(max.d, 6, 7)) != 0))){
    min.d <- min.d[-which(as.numeric(substr(min.d, 6, 7)) < 6)]
    max.d <- max.d[-which(as.numeric(substr(max.d, 6, 7)) > 6)]
  }
  
  range <- range(min.d, max.d) %>%
    substr(., 6, 10) %>%
    paste(., sep=" ", collapse=" to ") 
  
  return(range)
})

cowplot::plot_grid(plotlist = SST.plots, nrow = 1, labels = date.titles)

#my.plot <- ggplot() +
#  geom_jitter(fit.mpm@data, width = 0.2, pch = 16, alpha = 0.3, mapping = aes(x = iseason.ras, y = SST, colour = SST)) +
#  scale_colour_viridis_c(option = "B", limits = c(0, 35)) + 
#  theme(legend.position = "None") + 
#  xlab("Month") + 
#  ylab("SST (°C)") +
#  labs(title = "Monthly trends in sea surface temperature") 

############################################################################################################

import.river.data()
pretty.riv <- as(distance.any.ras, "SpatialPixelsDataFrame")
pretty.riv <- crop(pretty.riv, caspian.shape)

pretty.riv.pl <- ggplot() +
  gg(pretty.riv, aes(fill = (layer))) +
  gg(river.complex.lines, lwd = 2, aes(colour = Complex)) +
  gg(river.inlets, size = 3, aes(colour = Complex)) +
  scale_fill_gradientn(colours = spat.pal, name = "Distance from\nrivers (km)") +
  theme.standards +
  theme.spatials +
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 
############################################################################################################

contour.line <- aggregate(Bathymetry, fact=10) < -50 
contour.line[contour.line == 0] <- NA
contour.line <- rasterToPolygons(contour.line, dissolve = TRUE)
contour.line <- as(contour.line, "SpatialLines")
contour.line <- spTransform(contour.line, CRS(proj4string(fit.mpm)))
contour.line <- smoothr::smooth(contour.line, method = "ksmooth", smoothness = 5) #smoothen abit
contour.line <- crop(contour.line, caspian.shape)
contour <- as(contour, "SpatialPixelsDataFrame")
contour <- crop(contour, caspian.shape)

pretty.con.pl <- ggplot() +
  gg(contour, aes(fill = (layer))) +
  gg(contour.line, lwd = 2, colour = "grey20") +
  scale_fill_gradientn(colours = spat.pal, name = "Distance from\n50m contour (km)") +
  theme.standards +
  theme.spatials +
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 



############################################################################################################

#RESULTS AND FIGURES

length(unique(fit.mpm$id))
how.long.last <- function(data){
  data <- split(data, data$id) %>%
    lapply(., function(x){
      x$date <- as.Date(x$date, tryFormats = "%Y-%m-%d")
      return(as.numeric(max(x$date) - min(x$date)))
    })
  return(data)
}

round(mean(do.call(c, how.long.last(fit.mpm))), 2)
round(sd(do.call(c, how.long.last(fit.mpm))), 2)

#########################################################
wint.tags <- subset(fit.mpm, substr(datetagged, 4,5) %in% c("10", "11"))
sum.tags <- subset(fit.mpm, substr(datetagged, 4,5) %in% c("04", "05"))

round(mean(do.call(c, how.long.last(wint.tags))), 2)
round(sd(do.call(c, how.long.last(wint.tags))), 2)

round(mean(do.call(c, how.long.last(sum.tags))), 2)
round(sd(do.call(c, how.long.last(sum.tags))), 2)

round(nrow(wint.tags) / nrow(fit.mpm), 2)
round(nrow(sum.tags) / nrow(fit.mpm), 2)


table(fit.mpm$sex)[1] / nrow(fit.mpm)
table(fit.mpm$sex)[2] / nrow(fit.mpm)



#########################################################

cowplot::plot_grid(pretty.riv.pl, pretty.con.pl)


model.time$summary.fixed

model.time$summary.hyperpar
#range of spatial effect is low comapred to maximal extent fo caspian, therefore high ish clustering
abs(extent(caspian.shape)[1] - extent(caspian.shape)[2])
abs(extent(caspian.shape)[3] - extent(caspian.shape)[4])






#The total distance traveled by vessels varied widely between months, 
#ranging between 1004092 km in November to 1538251 km in April (Mean 1333218 SD 198896).
split(summary.output.paplan, summary.output.paplan$iseason) %>%
  lapply(., function(x){
    return(data.frame(iseason = unique(x$iseason), 
                      mid.date = unique(x$mid.date), 
                      total.dis = round(subset(x, zone == "total")$total.dis)))
  }) %>%
  do.call("rbind", .)
mean(subset(summary.output.paplan, zone == "total")$total.dis)
sd(subset(summary.output.paplan, zone == "total")$total.dis)
#There was less seasonal variability in the proportion of vessel travel which occured within each MPA plan (Figure XXX).
#Overall, unpenalized MPAs had a higher proportion of vessel travel occurring within their boundaries, and vessel density
#penalized MPAs had a lower proportion of vessel travel occurring within their boundaries.
#MPAs which were designed to cover 10% of the cumulative spatial ARS prediction contained, on average, 6.57% (SD 2.42) 
#in unpenalized plans and 2.16% (SD 0.74) in vessel penalized plans, of the total distance traveled by vessels within any given season.
#MPAs which were designed to cover 30% of the cumulative spatial ARS prediction contained, on average, 20.36% (SD 3.16)
#in unpenalized plans and 10.50% (SD 2.09 in vessel penalized plans, of the total distance traveled by vessels within any given season.
split(summary.output.paplan, summary.output.paplan$zone) %>%
  lapply(., function(x){
    return(data.frame(zone = unique(x$zone), mean = round(mean(x$perc.dis)*100, 2), sd = round(sd(x$perc.dis)*100, 2)))
  }) %>%
  do.call("rbind", .)
#Overall, the average travel speed of vessels did not vary widely between seasons (Mean 6.87 SD 0.35), 
#ranging between 6.03 kmh (SD 6.57) in October to 7.33 (SD 6.96) kmh in August
mean(subset(summary.output.paplan, zone == "total")$mean.speed.kmh)
sd(subset(summary.output.paplan, zone == "total")$mean.speed.kmh)
split(summary.output.paplan, summary.output.paplan$iseason) %>%
  lapply(., function(x){
    return(data.frame(iseason = unique(x$iseason), 
                      mid.date = unique(x$mid.date), 
                      mean = round(subset(x, zone == "total")$mean.speed.kmh, 2), 
                      sd = round(subset(x, zone == "total")$sd.speed, 2)))
  }) %>%
  do.call("rbind", .)
#However, the average travel speed of vessels within MPA plans did vary between seasons and were on average
#higher than overall vessel speeds.
#Vessel speeds (kmh) within MPAS were at their highest in July (Mean 12.91 SD 1.54) and their lowest in 
#October (Mean 6.16 SD 1.78)
split(summary.output.paplan, summary.output.paplan$iseason) %>%
  lapply(., function(x){
    x <- subset(x, zone != "total")
    return(data.frame(iseason = unique(x$iseason), 
                      mid.date = unique(x$mid.date), 
                      mean = round(mean(x$mean.speed.kmh), 2), 
                      sd = round(sd(x$mean.speed.kmh), 2)))
  }) %>%
  do.call("rbind", .)
#There was generally little difference between the temporal history of vessel speeds between each MPA categorization.
#MPAs which were designed to cover 10% of the cumulative spatial ARS prediction had a mean vessel travel travel speed of 
#9.54 kmh (SD 2.55) in unpenalized plans and 9.78 kmh (SD 3.54) in vessel penalized plans. MPAs which were designed 
#to cover 30% of the cumulative spatial ARS prediction had a mean vessel travel travel speed of 9.48 kmh (SD 1.31) in 
#unpenalized plans and 10.31 kmh (SD 2.30) in vessel penalized plans.





















covar.plot <- cowplot::plot_grid(pretty.riv.pl, pretty.con.pl)


coef.plot <- cowplot::plot_grid(fixed.plot, temporal.plot.full)

spatialmap <- cowplot::plot_grid(nrow = 4, ncol = max(fit.mpm$iseason.ras)+1, 
                                 plotlist = c(model.maps,
                                              heat.maps,
                                              priorit.maps.1,
                                              priorit.maps.2), byrow = T)

mpa.temp.plots <- cowplot::plot_grid(mpa.dis.perc.cover.plots, mpa.speed.plots, 
                                     cowplot::get_legend(mpa.dis.perc.cover.plots + theme(legend.position = "right")),
                                     rel_widths = c(1,1,0.3),
                                     nrow = 1, ncol = 3, byrow = T)

cowplot::save_plot("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Plots/Finalized/covar.plot.png",
                   plot = covar.plot, 
                   base_width = 8,
                   base_height = 6.5)

cowplot::save_plot("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Plots/Finalized/coef.plot.png",
                   plot = coef.plot, 
                   base_width = 8,
                   base_height = 6.5)

cowplot::save_plot("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Plots/Finalized/spatialmap.png",
                   plot = spatialmap, 
                   base_width = 12.5,
                   base_height = 8)

cowplot::save_plot("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Plots/Finalized/mpa.temp.plots.png",
                   plot = mpa.temp.plots, 
                   base_width = 8,
                   base_height = 6.5)





