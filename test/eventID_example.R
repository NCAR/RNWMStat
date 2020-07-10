# example script on how to run event identification and matching 
#
options(warn=-1)

rm(list=ls())

source("eventIdentification.R")
source("matchEvents.R")
source("plot_events_obs_model.R")

#file1 <- "02056900_obs_mod.csv"
#data1 <- read.csv(file1,header=TRUE)
#data1$X <- data1$PCP <- NULL
#data1$Date <- as.POSIXct(as.character(data1$Date),format="%Y-%m-%d %H:00:00",tz="UTC")

#for (gage0 in c("05056100","09306242","10172860","10291500","06186500","09386900","10243260")) {
#for (gage0 in c("01609000")) {
#for (gage0 in c("02310525")) {
#for (gage0 in c("04214500")) {
for (gage0 in c("02HC025")) {

print(paste0("------------- ", gage0, " -----------------"))

#dir0 <- "/glade/scratch/xfeng/NWMV3.0_Calib/objfun/retro_run/"
dir0 <- "/glade/scratch/xfeng/NWMV3.0_Calib/objfun/calib_run/"
data1 <- get(load(paste0(dir0,gage0,".Rdata")))
data1 <- as.data.frame(data1)
names(data1) <- c("site_no","Date","mod","obs")

# parameters for non-snow basins
nwinSpan1 <- 36
minEventDuration1 <- 12
maxRiseDuration1 <- 2*24
maxRecessionDuration1 <- 5*24
maxPeakDistCompound1 <- 5*24
maxDist1 <- 2*24

# parameters for seasonal snow basins
#nwinSpan1 <- 240
#minEventDuration1 <- 36
#maxRiseDuration1 <- 30*24
#maxRecessionDuration1 <- 60*24
#maxPeakDistCompound1 <- 15*24
#maxDist1 <- 7*24

# identify events for observed streamflow
print("identify events for observed streamflow")
threshPeak1 <- 0.98
threshFlowRange1 <- 0.3 
listObs <- eventIdentification(data1[,c("Date","obs")],nwinSpan=nwinSpan1,threshPeak=threshPeak1,threshFlowRange=threshFlowRange1, minEventDuration=minEventDuration1, maxRiseDuration=maxRiseDuration1,maxPeakDistCompound=maxPeakDistCompound1,maxRecessionDuration=maxRecessionDuration1)

# identify events for model streamflow
print("identify events for model streamflow")
threshPeak1 <- 0.90
threshFlowRange1 <- 0.3 
listMod <- eventIdentification(data1[,c("Date","mod")],nwinSpan=nwinSpan1,threshPeak=threshPeak1,threshFlowRange=threshFlowRange1,minEventDuration=minEventDuration1, maxRiseDuration=maxRiseDuration1,maxPeakDistCompound=maxPeakDistCompound1,maxRecessionDuration=maxRecessionDuration1)

# match observed events with model events
print("match observed events with model events")
eventsMatched <- matchEvents(data1[,c("Date","mod")],listMod[[1]], listMod[[2]], listObs[[2]], maxDist=maxDist1)

n1 <- sum(eventsMatched$match %in% c(1,2))
n2 <- sum(eventsMatched$match==3)
print(paste0("+++++ Number of events detected and matched: ", n1))
print(paste0("+++++ Number of events detected but missed by model: ", n2))

if (nrow(eventsMatched)>0) {
# peak timing error
time_err <- mean(abs(as.integer(difftime(eventsMatched$peak_mod,eventsMatched$peak_obs,units="hour"))))
print(paste0("timing error (hours): ", round(time_err,2)))

# peak bias (%)
obs_peak <- data1$obs[match(eventsMatched$peak_obs, data1$Date)]
mod_peak <- data1$mod[match(eventsMatched$peak_mod, data1$Date)]
peak_bias <- mean((mod_peak-obs_peak)/obs_peak*100)
print(paste0("peak bias (%): ", round(peak_bias,2)))

# volume bias (%)
ne <- nrow(eventsMatched)
volume_bias <- 0
for (i1 in 1:ne) {
  k1 <- match(eventsMatched$start_obs, data1$Date)
  k2 <- match(eventsMatched$end_obs, data1$Date)
  obs1 <- data1$obs[k1:k2]
  k1 <- match(eventsMatched$start_mod, data1$Date)
  k2 <- match(eventsMatched$end_mod, data1$Date)
  mod1 <- data1$mod[k1:k2]
  volume_bias <- volume_bias + (sum(mod1)-sum(obs1))/sum(obs1)*100
}
volume_bias <- volume_bias/ne
print(paste0("volume bias (%): ", round(volume_bias,2)))

# event-based objective function
obj1 <- 0.6*peak_bias + 0.4*volume_bias
print(paste0("event-based objective function: ", round(obj1,2)))

thresh1 <- quantile(data1$obs,0.9)
plotEventsObsMod(data1,listObs[[2]],listMod[[2]],gage0,thresh1)
}}
