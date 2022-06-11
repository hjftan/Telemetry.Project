library(sp)
library(raster)
library(ncdf4)
library(ggplot2)
library(tidyverse)
library(sf)
library("foieGras")
library(quantmod)
library(rgeos)

load("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Telemetry.clean.calculate.mpm.RData")

################################################################################################################################################
#this script 
#- uses output of script Telemetry.clean.calculate.mpm.R
#- creates and extracts environmental data
######################################################################################

my.proj <- "+proj=laea +lon_0=50.625 +lat_0=42.1915846 +datum=WGS84 +units=km +no_defs"

fit.mpm <- subset(fit.mpm, !is.na(x))
fit.mpm <- SpatialPointsDataFrame(coords = cbind(fit.mpm$x, fit.mpm$y), 
                                   data = fit.mpm[,-which(names(fit.mpm) %in% c("x", "y"))],
                                   proj4string = CRS("+init=EPSG:4326"))  #pretty sure its still 4326 here
fit.mpm <- spTransform(fit.mpm, CRS(my.proj))

##############
#define seasons which have bioloigcal relevance to seal ecology 

fit.mpm$season <- NA
fit.mpm$season[which(substr(fit.mpm$date, 6, 10) %in% 
                        substr(seq(from = as.Date("1066-01-01"), to = as.Date("1066-03-14"), by = 1), 6, 10))] <- 1 
#  1 = breeding
fit.mpm$season[which(substr(fit.mpm$date, 6, 10) %in% 
                        substr(seq(from = as.Date("1066-03-15"), to = as.Date("1066-04-16"), by = 1), 6, 10))] <- 2 
# 2 = moult# end moult at first looped tag

fit.mpm$season[which(substr(fit.mpm$date, 6, 10) %in% 
                        substr(seq(from = as.Date("1066-04-17"), to = as.Date("1066-06-07"), by = 1), 6, 10))] <- 3 
# 3 = early foraging # start foraging at end of moult, 4-17

fit.mpm$season[which(substr(fit.mpm$date, 6, 10) %in% 
                        substr(seq(from = as.Date("1066-06-08"), to = as.Date("1066-09-16"), by = 1), 6, 10))] <- 4 
# 4 = peak foraging 

fit.mpm$season[which(substr(fit.mpm$date, 6, 10) %in% 
                        substr(seq(from = as.Date("1066-09-17"), to = as.Date("1066-10-31"), by = 1), 6, 10))] <- 5
# 5 = late foraging 

fit.mpm$season[which(substr(fit.mpm$date, 6, 10) %in% 
                        substr(seq(from = as.Date("1066-11-01"), to = as.Date("1066-12-31"), by = 1), 6, 10))] <- 6 
# 6 = haul out before ice forms
fit.mpm$month <- as.numeric(substr(fit.mpm$date, 6, 7))
fit.mpm$year.month <- as.numeric(paste0(fit.mpm$year, substr(fit.mpm$date, 6, 7)))

#####################################################################################################################
#extract environmental data

SST <- c(list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2008", full.names = T, recursive = T),
         list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2009", full.names = T, recursive = T),
         list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2010", full.names = T, recursive = T),
         list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2011", full.names = T, recursive = T),
         list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2012", full.names = T, recursive = T),
         list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2013", full.names = T, recursive = T),
         list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2016", full.names = T, recursive = T),
         list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2017", full.names = T, recursive = T),
         list.files("F:/Project.Caspian.Hardrive/Anna.Data/SST.MUR.NEW/2018", full.names = T, recursive = T))


extract.sst.data <- function(data){
  
  data <- split(data, data$date)
  data <- lapply(data, function(x){
    
    my.date <- as.character(unique(str_replace_all(x$date, "-", "")))
    
    my.sst <- raster(SST[grepl(my.date, SST)])
    #my.sst <- mean(my.sst, na.rm = T), should only be 1
    
    my.crs <- CRS(proj4string(x))
    x <- spTransform(x, CRS(proj4string(my.sst)))
    x$SST <- raster::extract(my.sst, x)
    
    x <- spTransform(x, CRS(proj4string(my.crs)))
    
    return(x)
  }
  )
  data <- do.call("rbind", data)
  return(data)
}

fit.mpm <- extract.sst.data(fit.mpm)

################################################
#some NA for SST. these seem to arrise from a few records when seals were very close to land
#only a few entries so just estimate from the closest-time non na sst record for each indiviudal seal

length(which(is.na(fit.mpm$SST) == TRUE))
# NOT MANY SSTSO JUST set to SST when last correct, cant do disntance cause of seasonal variation
for ( i in which(is.na(fit.mpm$SST) == TRUE)){
  my.id <- fit.mpm$id[i]
  my.valid.data <- fit.mpm[which((fit.mpm$id == my.id) & (!is.na(fit.mpm$SST))),]
  fit.mpm$SST[i] <- my.valid.data$SST[min(which(abs(my.valid.data$idate - fit.mpm$idate[i]) ==
                                                  min(abs(my.valid.data$idate - fit.mpm$idate[i]))))]
  rm(my.id)
  rm(my.valid.data)
}
length(which(is.na(fit.mpm$SST) == TRUE))

########################################################
#extract river data

#The caspian sea is fed by many rivers. River mouths, especially ones which form river deltas, are known
#to be utilized by foraging seals. If Caspian seals utilize river mouths for foraging, we would expect 
#an estimated foraging metric to be negatively correlated with their distance to a given riven mouth.
# In addition, Caspian seals are known to eat migratory fish species 
#which migrate into or out of the Caspian to spawn, we may therefore expect their use of these 
#fresh water inflows to be seasonally dependent. 

#To map river inlets which feed the Caspian sea, we considered the 9 tributaries which provide the largest 
#amount of fresh water inflow, the Volga, Samur, Kura, Astarachay, Terek, Sulak, Ural, and rivers stemming from
#the Alborz mountain range.
#To calculate distance from each tributory, we calculated a distance raster from the rivers main water inlet.
# In cases where a main trbutory has many outlets, such as the Volga and Ural, We mapped each of the most distant outlets and generated
# a line feature between each of them, the subsequent distance raster was then calculate as a distance from
#the adjoining line. 
#Similarly, we grouped together river main tributarys when the tribuatory are in close proximity 
#and stem from same fresh water source, such as with the Terek and Sulak, Kura and Astarachay, Sefid-Rud and Haraz of the Alborz range.
# We also considered the Samur river source to be all river inlets ranging from the sSamur major river inlet
#to the Sabrancay inlet.

#In total we summarized these 8 main tributarys into 6 main tributarys

import.river.data <- function(){
  river.inlets <- data.frame(Complex = c("Volga", "Volga",
                                         "Terek-Sulak", "Terek-Sulak",
                                         "Samur-Sabrancay", "Samur-Sabrancay",
                                         "Kura-Astarachay", "Kura-Astarachay",
                                         "Safarud-Haraz", "Safarud-Haraz",
                                         "Ural", "Ural"),
                             River = c("Volga-1", "Volga-2",
                                       "Terek", "Sulak",
                                       "Samur", "Sabrancay",
                                       "Kura", "Astarachay",
                                       "Safarud", "Haraz",
                                       "Ural-1", "Ural-2"),
                             Lon = NA,
                             Lat = NA)
  
  #degree minute second format
  river.inlets$Lon[which(river.inlets$River == "Volga-1")] <- "47 30 17.76 E"
  river.inlets$Lat[which(river.inlets$River == "Volga-1")] <- "45 25 17.65 N"
  
  river.inlets$Lon[which(river.inlets$River == "Volga-2")] <- "49 53 59.48 E"
  river.inlets$Lat[which(river.inlets$River == "Volga-2")] <- "46 33 46.13 N"
  
  river.inlets$Lon[which(river.inlets$River == "Terek")] <- "47 33 36.40 E"
  river.inlets$Lat[which(river.inlets$River == "Terek")] <- "43 35 41.34 N"
  
  river.inlets$Lon[which(river.inlets$River == "Sulak")] <- "47 32 19.76 E"
  river.inlets$Lat[which(river.inlets$River == "Sulak")] <- "43 15 32.71 N"
  
  river.inlets$Lon[which(river.inlets$River == "Samur")] <- "48 29 11.89 E"
  river.inlets$Lat[which(river.inlets$River == "Samur")] <- "41 54 47.30 N"
  
  river.inlets$Lon[which(river.inlets$River == "Sabrancay")] <- "49 05 22.52 E"
  river.inlets$Lat[which(river.inlets$River == "Sabrancay")] <- "41 20 18.78 N"
  
  river.inlets$Lon[which(river.inlets$River == "Kura")] <- "49 23 45.17 E"
  river.inlets$Lat[which(river.inlets$River == "Kura")] <- "39 22 50.89 N"
  
  river.inlets$Lon[which(river.inlets$River == "Astarachay")] <- "48 52 51.70 E"
  river.inlets$Lat[which(river.inlets$River == "Astarachay")] <- "38 26 33.11 N"
  
  river.inlets$Lon[which(river.inlets$River == "Safarud")] <- "50 39 13.08 E"
  river.inlets$Lat[which(river.inlets$River == "Safarud")] <- "36 56 32.59 N"
  
  river.inlets$Lon[which(river.inlets$River == "Haraz")] <- "52 26 45.49 E"
  river.inlets$Lat[which(river.inlets$River == "Haraz")] <- "36 40 49.78 N"
  
  river.inlets$Lon[which(river.inlets$River == "Ural-1")] <- "50 48 41.76 E"
  river.inlets$Lat[which(river.inlets$River == "Ural-1")] <- "47 02 22.04 N"
  
  river.inlets$Lon[which(river.inlets$River == "Ural-2")] <- "52 46 14.94 E"
  river.inlets$Lat[which(river.inlets$River == "Ural-2")] <- "46 56 14.79 N"
  
  angle2dec <- function(angle) {
    angle <- as.character(angle)
    x <- do.call(rbind, strsplit(angle, split=' '))
    x <- apply(x, 1L, function(y) {
      y <- as.numeric(y)
      y[1] + y[2]/60 + y[3]/3600
    })
    return(x)
  }
  
  river.inlets$Lon <- angle2dec(river.inlets$Lon) 
  river.inlets$Lat <- angle2dec(river.inlets$Lat) 
  river.inlets <- SpatialPointsDataFrame(coords = cbind(river.inlets$Lon, river.inlets$Lat),
                                         data = river.inlets,
                                         proj4string = CRS("+init=EPSG:4326"))
  river.complex.lines <- list()
  for ( i in unique(river.inlets$Complex)){
    my.sub <- subset(river.inlets, Complex == i)
    my.sub <- spLines(as.matrix(my.sub[,c("Lon", "Lat")]@data), crs = CRS(proj4string(river.inlets)))
    my.sub <- SpatialLinesDataFrame(my.sub, data = data.frame(Complex = i))
    river.complex.lines[[i]] <- my.sub
  }
  river.complex.lines <- do.call("rbind", river.complex.lines)
  
  river.inlets <- spTransform(river.inlets, CRS(proj4string(fit.mpm)))
  river.complex.lines <- spTransform(river.complex.lines, CRS(proj4string(fit.mpm)))
  
  river.complex.lines$Complex.Short <- c("NW", "UW", "MW", "SW", "SE", "NE")
  river.inlets$Complex.Short <- rep(c("NW", "UW", "MW", "SW", "SE", "NE"), each = 2)
  
  assign("river.complex.lines", river.complex.lines, globalenv())
  assign("river.inlets", river.inlets, globalenv())
}
import.river.data()

distance.ras <- list()
for ( i in c("NW", "UW", "MW", "SW", "SE", "NE")){
  distance.ras[[i]] <- distanceFromPoints(raster(extent(caspian.shape) + 1000, resolution = 3, crs = CRS(proj4string(caspian.shape))), 
                                          spsample(subset(river.complex.lines, Complex.Short == i), 
                                                   gLength(subset(river.complex.lines, Complex.Short == i)), # so atleast 1 per km 
                                                   type = "regular"))
}

distance.ras <- stack(distance.ras)
distance.any.ras <- min(distance.ras)

extract.dis.inlet <- function(data){
  data$Dis.NW <- raster::extract(distance.ras[["NW"]], data)
  data$Dis.UW <- raster::extract(distance.ras[["UW"]], data)
  data$Dis.MW <- raster::extract(distance.ras[["MW"]], data)
  data$Dis.SW <- raster::extract(distance.ras[["SW"]], data)
  data$Dis.SE <- raster::extract(distance.ras[["SE"]], data)
  data$Dis.NE <- raster::extract(distance.ras[["NE"]], data)
  data$Dis.Any <- raster::extract(distance.any.ras, data)
  return(data)
}

fit.mpm <- extract.dis.inlet(fit.mpm)

#############################################################
#extract bathymetry and calculate and extract contour/isobath distances

Bathymetry <- raster("F:/Project.Caspian.Hardrive/Downloaded.Data/Bathymetry/gebco_2021_n47.85417_s36.15417_w45.90417_e54.875.tif")

Bathymetry <- Bathymetry + 27.711 # adjust for caspian sea level
Bathymetry[Bathymetry>0] <- 0

#200m isobath = contour https://royalsocietypublishing.org/doi/10.1098/rsos.210522#d1e2856
#bathmetric contours are potential important feeding areas for seals.
#we use -50m contour due to the caspian seals small size an typical diving capabilities

contour <- aggregate(Bathymetry, fact=10) < -50 # 
contour[contour == 0] <- NA
contour <- rasterToPolygons(contour, dissolve = TRUE)
contour <- as(contour, "SpatialLines")
contour <- spTransform(contour, CRS(my.proj))
contour <- spsample(contour, gLength(contour)*2, type = "regular") # 2 point per 1km
contour <- distanceFromPoints(raster(extent(caspian.shape) + 1000, 
                                      resolution = 1, 
                                      crs = CRS(proj4string(caspian.shape))),
                               contour)

slope <- terrain(Bathymetry, opt = "slope", unit = "degrees", neighbors = 4)

fit.mpm$bath <- raster::extract(Bathymetry, spTransform(fit.mpm, CRS(proj4string(Bathymetry))))
fit.mpm$slope <- raster::extract(slope, spTransform(fit.mpm, CRS(proj4string(slope))))
fit.mpm$bath.slope <- (1 + fit.mpm$bath + max(abs(fit.mpm$bath))) * fit.mpm$slope

fit.mpm$dis.con <- raster::extract(contour, spTransform(fit.mpm, CRS(proj4string(contour))))

#######################################################
#extract ice data

extract.ice.data <- function(data){
  
  data$ice.conc <- 999 # use 999 as missing ice value
  
  ice.data <- stack(list.files("F:/Project.Caspian.Hardrive/Downloaded.Data/asi.ice.telemetry.work/", full.names = T))
  
  my.crs <- CRS(proj4string(data))
  data <- spTransform(data, CRS(proj4string(ice.data)))
  
  my.dates <- unique(data$date)
  
  data <- split(data, data$date)
  
  for ( my.date in my.dates){
    
    my.ice <- ice.data[[which(grepl(gsub("-", "", my.date), names(ice.data)) == TRUE)]]
    
    if(length(my.ice) != 0){
      
      print(my.date)
      print(names(my.ice))
      
      data[[my.date]]$ice.conc <- raster::extract(my.ice, data[[my.date]])
      
    }
    
  }
  data <- do.call("rbind", data)
  data <- spTransform(data, my.crs)
  return(data)
  
}
fit.mpm <- extract.ice.data(fit.mpm)

hist(fit.mpm$ice.conc)
#255 and 999 are placeholders for land/nodata
fit.mpm$ice.conc[which(fit.mpm$ice.conc %in% c(255, 999))] <- 0
fit.mpm$ice.conc[which(is.na(fit.mpm$ice.conc))] <- 0
fit.mpm$ice.conc <- (fit.mpm$ice.conc - min(fit.mpm$ice.conc)) / (max(fit.mpm$ice.conc) - min(fit.mpm$ice.conc))

#####
#if on ice sometimes the movement of the ice/ amognst ice gives ARS, 
#correct for this by setting ARS/MPM to the average ARS value for each seal indepedently
ggplot() +
  geom_smooth(data = fit.mpm@data, aes(x=ice.conc, y = gmean.inv))

for( i in which(fit.mpm$ice.conc > 0.5)){
  fit.mpm$gmean[i] <- mean(subset(fit.mpm, index.id == fit.mpm$index.id[i])$gmean)
  fit.mpm$g.auc.stan[i] <- mean(subset(fit.mpm, index.id == fit.mpm$index.id[i])$g.auc.stan)
  fit.mpm$raw.mean[i] <- mean(subset(fit.mpm, index.id == fit.mpm$index.id[i])$raw.mean)
}

##################################################################################################################
#transform covariates and scale

scale.covariates <- function(data, scale.data){
  
  data$sq.dis.con <- sqrt(data$dis.con)
  
  data$sq.riv <- sqrt(data$Dis.Any)
  data$sq.riv.NW <- sqrt(data$Dis.NW)
  data$sq.riv.UW <- sqrt(data$Dis.UW)
  data$sq.riv.MW <- sqrt(data$Dis.MW)
  data$sq.riv.SW <- sqrt(data$Dis.SW)
  data$sq.riv.SE <- sqrt(data$Dis.SE)
  data$sq.riv.NE <- sqrt(data$Dis.NE)
  #data$sq.con <- sqrt(data$dis.contour.m)
  #data$sq.slope <- sqrt(data$slope)
  
  data$lg.riv <- log(data$Dis.Any)
  data$lg.riv.NW <- log(data$Dis.NW)
  data$lg.riv.UW <- log(data$Dis.UW)
  data$lg.riv.MW <- log(data$Dis.MW)
  data$lg.riv.SW <- log(data$Dis.SW)
  data$lg.riv.SE <- log(data$Dis.SE)
  data$lg.riv.NE <- log(data$Dis.NE)
  #data$lg.con <- log(data$dis.contour.m)
  #data$lg.slope <- log(data$slope)
  
  
  data$sq.dis.haul.any <- sqrt(data$dis.haul.any)
  data$sq.dis.haul.1 <- sqrt(data$dis.haul.1)
  data$sq.dis.haul.2 <- sqrt(data$dis.haul.2)
  data$sq.dis.haul.3 <- sqrt(data$dis.haul.3)
  data$sq.dis.haul.4 <- sqrt(data$dis.haul.4)
  data$sq.dis.haul.5 <- sqrt(data$dis.haul.5)
  data$sq.dis.haul.6 <- sqrt(data$dis.haul.6)
  data$sq.dis.haul.7 <- sqrt(data$dis.haul.7)
  data$sq.dis.haul.8 <- sqrt(data$dis.haul.8)
  data$sq.dis.haul.9 <- sqrt(data$dis.haul.9)
  
  data$cent.dis.con <- (data$dis.con - mean(c(scale.data$dis.con))) / sd(c(scale.data$dis.con))
  data$cent.sq.dis.con <- (data$sq.dis.con - mean(c(scale.data$sq.dis.con))) / sd(c(scale.data$sq.dis.con))
  
  data$cent.riv <- (data$Dis.Any - mean(c(scale.data$Dis.Any))) / sd(c(scale.data$Dis.Any))
  data$cent.sq.riv <- (data$sq.riv - mean(c(scale.data$sq.riv))) / sd(c(scale.data$sq.riv))
  data$cent.sq.riv.NW <- (data$sq.riv.NW - mean(c(scale.data$sq.riv.NW))) / sd(c(scale.data$sq.riv.NW))
  data$cent.sq.riv.UW <- (data$sq.riv.UW - mean(c(scale.data$sq.riv.UW))) / sd(c(scale.data$sq.riv.UW))
  data$cent.sq.riv.MW <- (data$sq.riv.MW - mean(c(scale.data$sq.riv.MW))) / sd(c(scale.data$sq.riv.MW))
  data$cent.sq.riv.SW <- (data$sq.riv.SW - mean(c(scale.data$sq.riv.SW))) / sd(c(scale.data$sq.riv.SW))
  data$cent.sq.riv.SE <- (data$sq.riv.SE - mean(c(scale.data$sq.riv.SE))) / sd(c(scale.data$sq.riv.SE))
  data$cent.sq.riv.NE <- (data$sq.riv.NE - mean(c(scale.data$sq.riv.NE))) / sd(c(scale.data$sq.riv.NE))
  
  data$cent.lg.riv <- (data$lg.riv - mean(c(scale.data$lg.riv))) / sd(c(scale.data$lg.riv))
  data$cent.lg.riv.NW <- (data$lg.riv.NW - mean(c(scale.data$lg.riv.NW))) / sd(c(scale.data$lg.riv.NW))
  data$cent.lg.riv.UW <- (data$lg.riv.UW - mean(c(scale.data$lg.riv.UW))) / sd(c(scale.data$lg.riv.UW))
  data$cent.lg.riv.MW <- (data$lg.riv.MW - mean(c(scale.data$lg.riv.MW))) / sd(c(scale.data$lg.riv.MW))
  data$cent.lg.riv.SW <- (data$lg.riv.SW - mean(c(scale.data$lg.riv.SW))) / sd(c(scale.data$lg.riv.SW))
  data$cent.lg.riv.SE <- (data$lg.riv.SE - mean(c(scale.data$lg.riv.SE))) / sd(c(scale.data$lg.riv.SE))
  data$cent.lg.riv.NE <- (data$lg.riv.NE - mean(c(scale.data$lg.riv.NE))) / sd(c(scale.data$lg.riv.NE))
  
  data$cent.sq.dis.haul.any <- (data$sq.dis.haul.any - mean(c(scale.data$sq.dis.haul.any))) / sd(c(scale.data$sq.dis.haul.any))
  data$cent.sq.dis.haul.1 <- (data$sq.dis.haul.1 - mean(c(scale.data$sq.dis.haul.1))) / sd(c(scale.data$sq.dis.haul.1))
  data$cent.sq.dis.haul.2 <- (data$sq.dis.haul.2 - mean(c(scale.data$sq.dis.haul.2))) / sd(c(scale.data$sq.dis.haul.2))
  data$cent.sq.dis.haul.3 <- (data$sq.dis.haul.3 - mean(c(scale.data$sq.dis.haul.3))) / sd(c(scale.data$sq.dis.haul.3))
  data$cent.sq.dis.haul.4 <- (data$sq.dis.haul.4 - mean(c(scale.data$sq.dis.haul.4))) / sd(c(scale.data$sq.dis.haul.4))
  data$cent.sq.dis.haul.5 <- (data$sq.dis.haul.5 - mean(c(scale.data$sq.dis.haul.5))) / sd(c(scale.data$sq.dis.haul.5))
  data$cent.sq.dis.haul.6 <- (data$sq.dis.haul.6 - mean(c(scale.data$sq.dis.haul.6))) / sd(c(scale.data$sq.dis.haul.6))
  data$cent.sq.dis.haul.7 <- (data$sq.dis.haul.7 - mean(c(scale.data$sq.dis.haul.7))) / sd(c(scale.data$sq.dis.haul.7))
  data$cent.sq.dis.haul.8 <- (data$sq.dis.haul.8 - mean(c(scale.data$sq.dis.haul.8))) / sd(c(scale.data$sq.dis.haul.8))
  data$cent.sq.dis.haul.9 <- (data$sq.dis.haul.9 - mean(c(scale.data$sq.dis.haul.9))) / sd(c(scale.data$sq.dis.haul.9))
  
  #data$cent.sq.con <- (data$sq.con - mean(c(scale.data$sq.con))) / sd(c(scale.data$sq.con))
  #data$cent.lg.con <- (data$lg.con - mean(c(scale.data$lg.con))) / sd(c(scale.data$lg.con))
  
  #data$cent.slope <- (data$slope - mean(c(scale.data$slope))) / sd(c(scale.data$slope))
  #data$cent.sq.slope <- (data$sq.slope - mean(c(scale.data$sq.slope))) / sd(c(scale.data$sq.slope))
  #data$cent.lg.slope <- (data$lg.slope - mean(c(scale.data$lg.slope))) / sd(c(scale.data$lg.slope))
  
  data$sq.slp <- sqrt(data$slope)
  data$cent.slp <- (data$slope - mean(c(scale.data$slope))) / sd(c(scale.data$slope))
  data$cent.sq.slp <- (data$sq.slp - mean(c(scale.data$sq.slp))) / sd(c(scale.data$sq.slp))
  data$cent.bath.slope <- (data$bath.slope - mean(c(scale.data$bath.slope))) / sd(c(scale.data$bath.slope))
  
  #data$cent.slp.mi10.bath <- (data$slp.mi10.bath - mean(c(scale.data$slp.mi10.bath))) / sd(c(scale.data$slp.mi10.bath))
  #data$cent.slp.mi20.bath <- (data$slp.mi20.bath - mean(c(scale.data$slp.mi20.bath))) / sd(c(scale.data$slp.mi20.bath))
  #data$cent.slp.mi30.bath <- (data$slp.mi30.bath - mean(c(scale.data$slp.mi30.bath))) / sd(c(scale.data$slp.mi30.bath))
  #data$cent.slp.mi40.bath <- (data$slp.mi40.bath - mean(c(scale.data$slp.mi40.bath))) / sd(c(scale.data$slp.mi40.bath))
  #data$cent.slp.mi50.bath <- (data$slp.mi50.bath - mean(c(scale.data$slp.mi50.bath))) / sd(c(scale.data$slp.mi50.bath))
  
  data$SST <- as.numeric(data$SST)
  data$cent.sst <- (data$SST - mean(c(scale.data$SST))) / sd(c(scale.data$SST))
  data$nrm.sst <- (data$SST - min(c(scale.data$SST))) / (max(c(scale.data$SST)) - min(c(scale.data$SST)))
  
  data$bath.abs <- abs(data$bath)
  data$bath.lg.abs <- log(data$bath.abs + 1) # plus 1 so above 0 
  
  data$cent.bath.abs <- (data$bath.abs - mean(c(scale.data$bath.abs))) / sd(c(scale.data$bath.abs))
  data$cent.bath.lg.abs <- (data$bath.lg.abs - mean(c(scale.data$bath.lg.abs))) / sd(c(scale.data$bath.lg.abs))
  
  return(data)
}
#data is data to be scaled, scale.data is data to scale by. useful to seprate for scaling prediciton data

fit.mpm <- scale.covariates(data = fit.mpm, scale.data = fit.mpm)


rm(list = ls()[-which(ls() %in% c("caspian.shape", "fit.mpm", "mesh",
                                  "slope.maxin.10",
                                  "slope.maxin.20",
                                  "slope.maxin.30",
                                  "slope.maxin.40",
                                  "slope.maxin.50",
                                  "distance.any.ras",
                                  "imma.haulouts",
                                  "dis.haulouts",
                                  "dis.haulouts.sep",
                                  "slope",
                                  "Bathymetry",
                                  "contour",
                                  "river.complex.lines", 
                                  "river.inlets"))])

#save.image("F:/Project.Caspian.Hardrive/RScripts/Telemetry.Work/Telemetry.create.extract.enviro.data.RData")



