library(sp)
library(raster)
library(ncdf4)
library(ggplot2)
library(tidyverse)
library(sf)
library(foieGras)
library(quantmod)

######################################################################################################
#this script 
#- imports raw satellite telemetry data from seal tags
#- reformats and cleans erroneous entries
#- fits state space model to standaridize track paths
#- fits movement persistance model to approximate area-restricted-search behaviours
#######################################################################################################
#load in raw data

list.files("Tag_Seals_details/csv.splits", 
           pattern = ".Individual.Details.csv",
           full.names = T)

Indiviudal.Details.2016 <- read.csv(list.files("Tag_Seals_details/csv.splits", 
                                               pattern = ".Individual.Details.csv",
                                               full.names = T)[1], stringsAsFactors = FALSE)
names(Indiviudal.Details.2016)[which(names(Indiviudal.Details.2016) == "Date")] <- "Date.Tagged"
Indiviudal.Details.2016$Seal.ID <- NULL

Indiviudal.Details.2017 <- read.csv(list.files("Tag_Seals_details/csv.splits", 
                                               pattern = ".Individual.Details.csv",
                                               full.names = T)[2], stringsAsFactors = FALSE)
names(Indiviudal.Details.2017)[which(names(Indiviudal.Details.2017) == "Date")] <- "Date.Tagged"
Indiviudal.Details.2017$Seal.ID <- NULL

Indiviudal.Details.2008.12 <- read.csv(list.files("Tag_Seals_details/csv.splits", 
                                                  pattern = ".Individual.Details.csv",
                                                  full.names = T)[3], stringsAsFactors = FALSE)
names(Indiviudal.Details.2008.12)[which(names(Indiviudal.Details.2008.12) == "Tagging.date")] <- "Date.Tagged"
names(Indiviudal.Details.2008.12)[which(names(Indiviudal.Details.2008.12) == "Wt..kg.")] <- "Weight..Kg."
names(Indiviudal.Details.2008.12)[which(names(Indiviudal.Details.2008.12) == "Body.length..cm.")] <- "Length..cm."
names(Indiviudal.Details.2008.12)[which(names(Indiviudal.Details.2008.12) == "Tag.no")] <- "Tag.ID"
Indiviudal.Details.2008.12$Tag.ID <- gsub("[^[:digit:]]", "", Indiviudal.Details.2008.12$Tag.ID)
Indiviudal.Details.2008.12$BMI. <- NULL
Indiviudal.Details.2008.12$Capture.Location <- NA

Indiviudal.Details.All <- rbind(Indiviudal.Details.2008.12, Indiviudal.Details.2016, Indiviudal.Details.2017)
Indiviudal.Details.All$Weight..Kg.
Indiviudal.Details.All$Length..cm.[which(Indiviudal.Details.All$Length..cm. < 10)] <- 
  Indiviudal.Details.All$Length..cm.[which(Indiviudal.Details.All$Length..cm. < 10)] * 100
Indiviudal.Details.All$Girth..cm.[which(Indiviudal.Details.All$Girth..cm. < 10)] <- 
  Indiviudal.Details.All$Girth..cm.[which(Indiviudal.Details.All$Girth..cm. < 10)] * 100
#some recorded in meters

###########################################################################################
#load filtered raw ssm data 

list.files("Filtered_raw_data",
           pattern = ".csv", full.names = T)

Filtered.Raw.SSM.2008 <- read.csv(list.files("Filtered_raw_data",
                                             pattern = ".csv", full.names = T)[1], stringsAsFactors = FALSE)
Filtered.Raw.SSM.2009 <- read.csv(list.files("Filtered_raw_data",
                                             pattern = ".csv", full.names = T)[2], stringsAsFactors = FALSE)
Filtered.Raw.SSM.2010 <- read.csv(list.files("Filtered_raw_data",
                                             pattern = ".csv", full.names = T)[3], stringsAsFactors = FALSE)
Filtered.Raw.SSM.2011 <- read.csv(list.files("Filtered_raw_data",
                                             pattern = ".csv", full.names = T)[4], stringsAsFactors = FALSE)
Filtered.Raw.SSM.2012 <- read.csv(list.files("Filtered_raw_data",
                                             pattern = ".csv", full.names = T)[5], stringsAsFactors = FALSE)
Filtered.Raw.SSM.2016.1 <- read.csv(list.files("Filtered_raw_data",
                                               pattern = ".csv", full.names = T)[6], stringsAsFactors = FALSE)
Filtered.Raw.SSM.2016.2 <- read.csv(list.files("Filtered_raw_data",
                                               pattern = ".csv", full.names = T)[7], stringsAsFactors = FALSE)
Filtered.Raw.SSM.2017 <- read.csv(list.files("Filtered_raw_data",
                                             pattern = ".csv", full.names = T)[8], stringsAsFactors = FALSE)
#the pos number is just an ordered factor for the date sequence, 2017 doesnt have it for some reason
Filtered.Raw.SSM.2008$Pos.number <- NULL
Filtered.Raw.SSM.2009$Pos.number <- NULL
Filtered.Raw.SSM.2010$Pos.number <- NULL
Filtered.Raw.SSM.2011$Pos.number <- NULL
Filtered.Raw.SSM.2012$Pos.number <- NULL
Filtered.Raw.SSM.2016.1$Pos.number <- NULL
Filtered.Raw.SSM.2016.2$Pos.number <- NULL

Filtered.Raw.SSM.All <- rbind(Filtered.Raw.SSM.2008, Filtered.Raw.SSM.2009, Filtered.Raw.SSM.2010, Filtered.Raw.SSM.2011,
                              Filtered.Raw.SSM.2012, Filtered.Raw.SSM.2016.1, Filtered.Raw.SSM.2016.2, Filtered.Raw.SSM.2017)

Filtered.Raw.SSM.All.WTagDets <- merge(Filtered.Raw.SSM.All, Indiviudal.Details.All, by.x = "tag.id", by.y = "Tag.ID")
Filtered.Raw.SSM.All.WTagDets <- SpatialPointsDataFrame(coords = cbind(Filtered.Raw.SSM.All.WTagDets$Lon, Filtered.Raw.SSM.All.WTagDets$Lat),
                                                        data = Filtered.Raw.SSM.All.WTagDets[,-which(names(Filtered.Raw.SSM.All.WTagDets) %in% c("Lon", "Lat"))])

data.ssm <- Filtered.Raw.SSM.All.WTagDets
data.ssm <- data.frame(id = data.ssm$tag.id, 
                       date = as.character(data.ssm$Date), 
                       lc = data.ssm$LC, lon = data.ssm@coords[,1], 
                       lat = data.ssm@coords[,2])

data.ssm$date <- as.POSIXlt(paste(data.ssm$date, data.ssm$time), format = "%d/%m/%Y %H:%M", tz = "GMT")

data.ssm[which(data.ssm$lc == -1),]$lc <- "A"
data.ssm[which(data.ssm$lc == -2),]$lc <- "B"
data.ssm[which(data.ssm$lc == 4),]$lc <- 3

number.pings.recorded <- sapply(split(data.ssm, data.ssm$id), function(x){return(nrow(x))})

data.ssm <- data.ssm[-which(data.ssm$id %in% names(number.pings.recorded[which(number.pings.recorded < 10)])),]

#################################################################################################
#fix duplicate entries by spreading them across an even/fixed time interval within each day
#duplicate entries seem to arise from when daily summaries were returned from the telelmetry data
#wont cause issues about psuedo-like replication because daily averages are calculated later

fix.duplicate.entries <- function(data.list){
  
  data.list <- lapply(data.list, function(x){
    
    dup.dates <- names(which(table(as.character(x$date)) != 1))
    
    for ( i in dup.dates) {
      
      dup.indexs <- which(as.character(x$date) == as.character(i))
      print(x[dup.indexs,])
      
      newdates <- seq(from = as.POSIXlt(paste(unique(x$date[dup.indexs]), "00:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT"),
                      to = as.POSIXlt(paste(unique(x$date[dup.indexs]), "23:59:59"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT"),
                      length.out = length(x$date[dup.indexs]) + 2) # plus two so take middle pieces
      newdates <- newdates[c(-1, -length(newdates))]
      
      x$date[dup.indexs] <- newdates
    
    }
    return(x)
  }
  ) # fix duplicate entries
  data.list <- do.call("rbind", data.list)
  row.names(data.list) <- NULL
  
  return(data.list)
}
data.ssm <- 
  split(data.ssm, data.ssm$id) %>%
  fix.duplicate.entries(.)

#################################################################################################

check.gaps <- split(data.clean, data.clean$id)

check.gaps <- lapply(check.gaps, function(x){
  my.dates <- sort(x$date)
  gap <- my.dates[2] - my.dates[1]
  
  for( i in 2:(length(my.dates)-1)){
    gap <- c(gap, my.dates[i + 1] - my.dates[i])
  }
  return(gap)
  
})
check.gaps.m <- lapply(check.gaps, function(x){
  return(((as.numeric(max(x))/60)/60)/24)})

check.gaps.m[which(check.gaps.m > 21)]

data.clean[which(data.clean$id %in% names(check.gaps.m[which(check.gaps.m > 21)])),]
data.clean$id[which((data.clean$id == 57039) & (data.clean$date < "2010-12-31"))] <- 57038
# 57039 is a tagid that got reused. replaced with 57038
replacement <- subset(Indiviudal.Details.All, Tag.ID == 57039)
replacement$Tag.ID <- 57038
Indiviudal.Details.All <- rbind(Indiviudal.Details.All, replacement)
rm(replacement)

data.clean <- data.clean[-which(data.clean$id %in% names(check.gaps.m[which(check.gaps.m > 21)])),]

####################################################################################################
#fit ssm to standardize the tracks
fit.ssm <- fit_ssm(data.clean, 
               model = "rw", # crw too many artefacts, and valid lots of smal lscale movmeents so use rw
               time.step = 6,
               #6 hour timestep so i can calculate daily averages of behaviour without losing small 
               #variations in movement. Important so smaller scale variations in space use are not lost
               #when i calculate the behavioural metric
               vmax = 8,
               min.dt = 30, # most dives are about a min for caspian seals, so need min dt to be less than that
               ang = c(50, 90),
               distlim = c(2500, 7500),
               control = ssm_control(verbose = 1, se = FALSE))

fit.mpm <- fit_mpm(fit.ssm, what = "predicted", model = "mpm", control = mpm_control(verbose = 0))
#calculate the movement persistance behavourial metric 

####################################################################################################

fitted.vals <- data.frame(grab(fmp, what = "fitted", as_sf = FALSE))
data.vals <- data.frame(grab(fmp, what = "data", as_sf = FALSE))
rownames(fitted.vals) <- NULL
rownames(data.vals) <- NULL

fitted.vals <- split(fitted.vals, fitted.vals$id)
data.vals <- split(data.vals, data.vals$id)
merge.vals <- list()
for(i in 1:length(data.vals)){
  merge.vals[[i]] <- merge(data.vals[[i]], fitted.vals[[i]], by = "date")
  
}
merge.vals <- do.call("rbind", merge.vals)
unique(merge.vals$id.x == merge.vals$id.y)
merge.vals <- merge.vals[,-which(names(merge.vals) %in% c("idtid", "tid", "dt", "id.y"))]
names(merge.vals)[which(names(merge.vals) == "id.x")] <- "id"

fit.mpm <- merge.vals

############################################################################################
#these tracks are too short i.e the seal didnt move very far or didnt last long enough
#to be a reliable record for foraging behaviour
bad.tags <- c(90776, 79888, 159405, 159385, 159383, 159381, 159380)

#these tags have erroneous looking movement persistance 
#estimates, clearly some problem in the calculation of mpm
#issue can be fixed for some by altering the base parameters for the ssm
#however, this  leads to other tags returning erroneous looking movement persistance
#this seems like the fewest I can remove whilst keeping tags which lasted a longer period of time
weird.tags <- c(57048, 57047, 57046, 57045, 57044, 57043, 57042, 57041, 57040, 57039, 57038, 57032, 57030,57026)
weird.tags.fmp <- subset(fit.mpm, id %in% as.character(weird.tags))

ggplot() +
  geom_line(data = subset(weird.tags.fmp, id == 57026), aes( x= date, y= g))
ggplot() +
  geom_line(data = subset(weird.tags.fmp, id == 57047), aes( x= date, y= g))

#nremove bad  tags
fit.mpm <- fit.mpm[-which(fit.mpm$id %in% bad.tags),]
fit.mpm <- fit.mpm[-which(fit.mpm$id %in% weird.tags),]

################################################################################################
#create indexs for modelling. e.g. sex, individual, data, year etc..

fit.mpm$indexid <- NA
for ( i in 1:length(unique(fit.mpm$id))) {
  fit.mpm$indexid[which(fit.mpm$id == unique(fit.mpm$id)[i])] <- i 
}

fit.mpm$indexdate <- NA
for ( i in 1:length(unique(fit.mpm$date))) {
  fit.mpm$indexdate[which(fit.mpm$date == sort(unique(fit.mpm$date)[i]))] <- i 
}

fit.mpm$year <- substr(fit.mpm$date, 1,4)
fit.mpm$iyear <- NA
for (i in 1:length(sort(unique(fit.mpm$year)))) {
  fit.mpm$iyear[which(fit.mpm$year == sort(unique(fit.mpm$year))[i])] <- i
}

fit.mpm$idatetagged <- NA
for ( i in 1:nrow(date.match)) {
  fit.mpm$idatetagged[which(paste0(substr(fit.mpm$datetagged, 1, 2), "/", substr(fit.mpm$datetagged, 4, 5)) == #incosnsitent format . or /
                              date.match[i,2])] <- as.numeric(date.match[i,1])
}

fit.mpm$indexsex <- NA
fit.mpm$indexsex[which(fit.mpm$sex == "M")] <- 1
fit.mpm$indexsex[which(fit.mpm$sex == "F")] <- 2

###################################################################################################

fit.mpm <- SpatialPointsDataFrame(coords = data.frame(cbind(fit.mpm$x, fit.mpm$y)), data = fit.mpm)

#we want to have a model effect which represents a annual "loop" i.e. tethered start and end points, 
# to do that its easier to model if each individual is only a member of a single loop year.
#this means the loop cant be 1st jan to 31 dec as a few datapoints cross over the christmas period.
#so we must define a doifferent start end point

check.df <- data.frame(start = substr(as.character(date.check.df$min.date), 6, 10),
                       end = substr(as.character(date.check.df$max.date), 6, 10))
ggplot() +
  geom_point(data = check.df, aes(x = start, y = end))
#earleist start is 3.26 but its only a few data points
subset(fit.mpm, substr(fit.mpm$date, 6, 10) == "03-26")
#i think better to start on the 04-17 
subset(fit.mpm, substr(fit.mpm$date, 6, 10) == "04-17")

test.date.seq <- seq(from = as.Date("1999-04-17"), to = as.Date("2000-04-16"), by = 1)
test.date.seq <- substr(test.date.seq, 6, 10)
ggplot() +
  geom_line(data = fit.mpm@data, aes(x=substr(fit.mpm$date, 6, 10), y = g, group = id)) +
  scale_x_discrete(limits = test.date.seq) +
  theme(axis.text.x = element_text(angle = 90))
#starting on 17th means only 1 point crosses over the anual loop
start.of.year <- "04-17"

early.starts <- unique(subset(fit.mpm, substr(fit.mpm$date, 6, 10) == "03-26")$id)
early.starts <- subset(fit.mpm, id %in% early.starts)
early.starts <- split(early.starts, early.starts$id)
early.starts <- lapply(early.starts, function(x){
  df <- data.frame(id = unique(x$id), start = min(x$date), end = max(x$date))
  return(df)
})
#the 1 ooint is id 57039, which only lastest a couple weeks. remove so doesnt cause problems
fit.mpm <- subset(fit.mpm, id != 57039)

############################################################################################
#create daily summaries of the mpm metric and fix idate = 1 to 17th of april
#summaries include e.g. daily mean and AUC-like approach
create.daily.averages.mpm <- function(data.mpm){
  
  my.data.simp <- data.frame(id = rep(unique(data.mpm$id), each = (365 / MODEL.INTERVAL)),
                             idate = rep(seq(from = 1, to = (365 / MODEL.INTERVAL)), times = length(unique(data.mpm$id))))
  my.data.simp$year <- NA # REMMEEBR A FAIR FEw data cross over december. will be a problem
  #if i wanted to fit year as an effect
  my.data.simp$iyear <- NA
  my.data.simp$date <- NA
  my.data.simp$x <- NA
  my.data.simp$y <- NA
  my.data.simp$mean <- NA
  my.data.simp$raw.mean <- NA
  my.data.simp$g.auc <- NA
  my.data.simp$g.auc.stan <- NA
  
  my.data.simp.list <- list()
  
  for ( i in unique(my.data.simp$id)) {
    print(i)
    
    my.dates <- subset(data.mpm, id == i)$date
    if (sum(my.dates < as.Date(paste0(max(as.numeric(substr(my.dates, 1, 4))), "-", start.of.year))) > 0){
      my.data.simp$year[which(my.data.simp$id == i)] <- max(as.numeric(substr(my.dates, 1, 4))) - 1
    } 
    if (sum(my.dates > as.Date(paste0(max(as.numeric(substr(my.dates, 1, 4))), "-", start.of.year))) > 0){
      my.data.simp$year[which(my.data.simp$id == i)] <- max(as.numeric(substr(my.dates, 1, 4)))
    } # no cross overs the start of year so wont be any double assignment
    if((sum(my.dates > as.Date(paste0(max(as.numeric(substr(my.dates, 1, 4))), "-", start.of.year))) > 0) &
       (sum(my.dates < as.Date(paste0(max(as.numeric(substr(my.dates, 1, 4))), "-", start.of.year))) > 0)){
      print("ERROR ERROR DOUBLRE ASIGNMENT, THIS ID LOOPS OVER START OF YEAR STOP AND CHECK CODE")
    }
    
    date.seq <- seq(from = as.POSIXct(paste0(min(as.numeric(unique(subset(data.mpm, id == i)$year))), "-", start.of.year), format = "%Y-%m-%d", tz = "GMT"),
                    to = as.POSIXct(paste0(min(as.numeric(unique(subset(data.mpm, id == i)$year))) + 1, "-", start.of.year), format = "%Y-%m-%d", tz = "GMT") - 1,
                    by = (((60*60)*24) * MODEL.INTERVAL))
    
    if(length(which(substr(date.seq, 6, 10) == "02-29")) != 0){
      date.seq <- date.seq[-which(substr(date.seq, 6, 10) == "02-29")]
      # if leap year # just delete it
    }
    
    my.data.simp$date[which(my.data.simp$id == i)] <- as.character(date.seq)
    
    my.sub <- cbind(data.mpm[which(data.mpm$id == i),]@data, data.mpm[which(data.mpm$id == i),]@coords)
    
    my.sub$short.date <- substr(my.sub$date, 1, 10)
    
    #for UTC
    time.df <- data.frame(times = seq(from = as.POSIXlt(paste(min(substr(my.sub$date,1,10)), "00:00:00"), tz = "UTC"),
                                      to = as.POSIXlt(paste(max(substr(my.sub$date,1,10)), "23:59:59"), tz = "UTC"), by = 60))
    
    time.df <- merge(x = time.df, y = my.sub[,which(names(my.sub) %in% c("date", "g", "x", "y"))], by.x = "times", by.y = "date", 
                     all.x = TRUE)
    
    #start and end of each day always set to first value recorded or last, otherwise the inteprolate dpoesnt work
    
    if(is.na(time.df$g[1])) {  time.df$g[1] <- time.df$g[which(!is.na(time.df$g))[1]] }
    if(is.na(time.df$g[nrow(time.df)])) {  time.df$g[nrow(time.df)] <- time.df$g[which(!is.na(time.df$g))[length(which(!is.na(time.df$g)))]] }
    
    if(is.na(time.df$x[1])) {  time.df$x[1] <- time.df$x[which(!is.na(time.df$x))[1]] }
    if(is.na(time.df$x[nrow(time.df)])) {  time.df$x[nrow(time.df)] <- time.df$x[which(!is.na(time.df$x))[length(which(!is.na(time.df$x)))]] }
    
    if(is.na(time.df$y[1])) {  time.df$y[1] <- time.df$y[which(!is.na(time.df$y))[1]] }
    if(is.na(time.df$y[nrow(time.df)])) {  time.df$y[nrow(time.df)] <- time.df$y[which(!is.na(time.df$y))[length(which(!is.na(time.df$y)))]] }
    
    #  plot(time.df$times, time.df$g)
    
    time.df$g <- zoo::na.approx(time.df$g)
    
    time.df$x <- zoo::na.approx(time.df$x)
    
    time.df$y <- zoo::na.approx(time.df$y)
    
    time.df$date <- substr(time.df$times, 1, 10)
    
    my.sub.id <- my.data.simp[which(my.data.simp$id == i),]
    my.sub.id <- split(my.sub.id, my.sub.id$idate)
    
    my.sub.id <- lapply(my.sub.id, FUN = function(j){
      
      if(MODEL.INTERVAL == 1) {
        
        index.my.sub <- which(my.sub$short.date == j$date)
        
        index.my.time <- which(time.df$date == j$date)
        
      } else{
        print("warning, if i want model interval not to be 1, runs slower and need the asDate, see this chunk")
        index.my.sub <- which((my.sub$short.date >= j$date) & (my.sub$short.date <= (as.Date(j$date) + (MODEL.INTERVAL - 1))))
        
        index.my.time <- which((time.df$date >= j$date) & 
                                 (time.df$date <= (as.Date(j$date) + (MODEL.INTERVAL - 1))))
      }
      
      if ((length(index.my.sub) != 0) | (length(index.my.time) != 0)) {
        
        j$x <- mean(time.df$x[index.my.time], na.rm = T) 
        j$y <- mean(time.df$y[index.my.time], na.rm = T) 
        
        #raw valjust from real data, same as old
        j$raw.mean <- mean(my.sub$g[index.my.sub], na.rm = T) 
        
        #mean from sampled data
        j$mean <- mean(time.df$g[index.my.time], na.rm = T) 
        
        #area under curve approach.
        #higher equals more time foraging, lower equals less time foraging
        #and stadnardized to max foraging for a model itnerval. e.g. value 1 for 24 hours
        j$g.auc <- MESS::auc(1:length(index.my.time), time.df$g[index.my.time], type = "linear") 
        j$g.auc.stan <- j$g.auc / ((24*60) * MODEL.INTERVAL)
        
      }
      return(j)
    }
    )
    
    my.data.simp.list[[i]] <- do.call("rbind", my.sub.id)
    
  }
  
  my.data.simp <- do.call("rbind", my.data.simp.list)
  return(my.data.simp)
  row.names(my.data.simp) <- NULL
  
}

fit.mpm <- create.daily.averages.mpm(fit.mpm)
#################################################################################################

names(Indiviudal.Details.All) <- c("id", "datetagged", "sex", "weightkg", "lengthcm", "girthcm", "capturelocation", "tagdateid")
fit.mpm <- merge(fit.mpm, Indiviudal.Details.All, by.x = "id", by.y = "id")
names(fit.mpm)[which(names(fit.mpm) == "mean")] <- "gmean" # otherewise name gets lost in prediciton

fit.mpm$index.id <- NA
for (i in 1:length(unique(fit.mpm$id))){
  fit.mpm$index.id[which(fit.mpm$id == unique(fit.mpm$id)[i])] <- i
}

#spo rmemebr our "year" starts in april
head(cbind(fit.mpm$idate, fit.mpm$date))
tail(cbind(fit.mpm$idate, fit.mpm$date))
date.match <- cbind(fit.mpm$idate, paste0(substr(fit.mpm$date, 9, 10), "/", substr(fit.mpm$date, 6, 7)))
date.match <- unique(date.match)
colnames(date.match) <- c("idate", "tagdate")

#save.image("Telemetry.clean.calculate.mpm.RData")


