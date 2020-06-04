rm(list=ls())

library(data.table)

source("eventIdentification.R")

# list of test basins
dfBasinFinal <- get(load("data/final_list_basins_attributes_short.Rdata"))
dfBasinFinal$snow <- 0
ix1 <- which(dfBasinFinal$gage %in% c("06289000","09081600","09492400")) 
dfBasinFinal$snow[ix1] <- 1

#for (i in 1:nrow(dfBasinFinal)) {
for (i in 1) {
gage <- dfBasinFinal$gage[i]

# model simulation 
file1 <- paste0("data/model/",gage,"_obs_mod.csv")
if (!file.exists(file1)) next
data1 <- read.csv(file1,header=TRUE)
data1$X <- data1$obs <- NULL
data1$Date <- as.POSIXct(as.character(data1$Date),format="%Y-%m-%d %H:00:00",tz="UTC")
data1 <- as.data.table(data1)

if (dfBasinFinal$snow[i]==0) {
nwinSpan1 <- 36
minEventDuration1 <- 12
maxRiseDuration1 <- 2*24
maxRecessionDuration1 <- 5*24
} else {
nwinSpan1 <- 168
minEventDuration1 <- 36
maxRiseDuration1 <- 30*24
maxRecessionDuration1 <- 60*24
}

list1 <- eventIdentification(data=data1,nwinSpan=nwinSpan1,
  minEventDuration=minEventDuration1, maxRiseDuration=maxRiseDuration1,
  maxRecessionDuration=maxRecessionDuration1)

events <- list1[[1]]
events[,t1:=as.integer(difftime(peak,start,units="hours"))]
events[,t2:=as.integer(difftime(end,peak,units="hours"))]
print(paste0(gage," ",nrow(events)))
print(range(events$t1))
print(range(events$t2))
}
