rm(list=ls())

workDir <- "/glade/work/ariz01/nwm3/runoff/"

source(paste0(workDir,"eventIdentification.R"))

# streamflow data
obsTestbasins <- get(load(paste0(workDir,"data/flow_obs_testbasins.Rdata")))

# list of test basins
dfBasinFinal <- get(load(paste0("data/final_list_basins_attributes_short.Rdata")))
dfBasinFinal$snow <- 0
ix1 <- which(dfBasinFinal$gage %in% c("06289000","09081600","09492400")) 
dfBasinFinal$snow[ix1] <- 1

for (i in 1:nrow(dfBasinFinal)) {
#for (i in 4) {
gage <- dfBasinFinal$gage[i]

# observation
obs <- subset(obsTestbasins, site_no == gage)
obs$site_no <- NULL

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

events <- eventIdentification(data=obs,tag=paste0(gage,"_obs"),nwinSpan=nwinSpan1,
  minEventDuration=minEventDuration1, maxRiseDuration=maxRiseDuration1,
  maxRecessionDuration=maxRecessionDuration1)

events[,t1:=as.integer(difftime(peak,start,units="hours"))]
events[,t2:=as.integer(difftime(end,peak,units="hours"))]
print(paste0(gage," ",nrow(events)))
#print(range(events$t1))
#print(range(events$t2))
}
