rm(list=ls())

library(ggplot2)
library(data.table)

workDir <- "/glade/work/ariz01/nwm3/runoff/"
prob_peak <- 0.9

source(paste0(workDir,"eventIdentification.R"))

# streamflow data
obsTestbasins <- get(load(paste0(workDir,"data/flow_obs_testbasins.Rdata")))

# list of test basins
dfBasinFinal <- get(load(paste0("data/final_list_basins_attributes_short.Rdata")))
dfBasinFinal$snow <- 0
ix1 <- which(dfBasinFinal$gage %in% c("06289000","09081600","09492400"))
dfBasinFinal$snow[ix1] <- 1

for (i in 1:nrow(dfBasinFinal)) {
#for (i in 12) {
#if (i==4) next
gage <- dfBasinFinal$gage[i]
#gage <- "05525500"

# observation
data <- subset(obsTestbasins, site_no == gage)
data$site_no <- NULL
names(data) <- c("time","value")
thresh1 <- quantile(data$value, prob_peak)

# snow and non-snow basins
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

# slow, groundwater basins in FL
if (gage %in% c("02310947","02479300")) {
  maxRiseDuration1 <- 30*24
  maxRecessionDuration1 <- 30*24
}

# relatively slow basins in IL
if (gage %in% c("05584500","05525500")) {
  maxRiseDuration1 <- 5*24
  maxRecessionDuration1 <- 8*24
}

list1 <- eventIdentification(data,nwinSpan=nwinSpan1, probPeak=prob_peak,
  minEventDuration=minEventDuration1, maxRiseDuration=maxRiseDuration1,
  maxRecessionDuration=maxRecessionDuration1)
eventsAll <- list1[[1]]
dataAll <- list1[[2]]

eventsAll[,t1:=as.integer(difftime(peak,start,units="hours"))]
eventsAll[,t2:=as.integer(difftime(end,peak,units="hours"))]
print(paste0(gage," ",nrow(eventsAll)))

outfile <- paste0("events/events_",gage,"_obs.Rdata")
if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile),recursive=TRUE)
save(eventsAll,file=outfile)

#plot hydrographs and peaks (in half-year chunks)
plot1 <- TRUE
if (plot1) {
years <- as.integer(unique(format(dataAll$time,"%Y")))
kk <- 0
for (y1 in years) {
for (k1 in 1:2) {
  if (k1==1) {
    dates <- seq(as.POSIXct(paste0(y1-1,"1001"),format="%Y%m%d"),
                 as.POSIXct(paste0(y1,"0401"),format="%Y%m%d"),by="hour")
  } else {
    dates <- seq(as.POSIXct(paste0(y1,"0401"),format="%Y%m%d"),
                 as.POSIXct(paste0(y1,"1001"),format="%Y%m%d"),by="hour")
  }
  data2 <- subset(dataAll, time %in% dates)
  #if (sum(peaksAll$time %in% data2$time)==0) next
  if (sum(eventsAll$peak %in% data2$time)==0) next

  # fill in missing times with NaN
  dates1 <- dates[!dates %in% data2$time]
  if (length(dates1)>0) {
  data2 <- rbind(data2,data.frame(time=dates1, value=NA, smooth=NA))
  data2 <- data2[order(time),]
  }

  peaks <- subset(data2, time %in% eventsAll$peak)
  starts <- subset(data2, time %in% eventsAll$start)
  ends <- subset(data2, time %in% eventsAll$end)

  kk <- kk + 1
  gg1 <- ggplot(data=data2,aes(time,value)) +
      geom_line(color="darkgrey") +
      geom_line(aes(time,smooth),color="black",size=0.4)+
      geom_point(data=peaks,aes(time,value),color="red",shape=2,size=1) +
      geom_point(data=starts,aes(time,value),color="blue",shape=8,size=1) +
      geom_point(data=ends,aes(time,value),color="purple",shape=1,size=1) +
      geom_hline(yintercept=thresh1,linetype="dashed",color="black",size=0.5) +
      theme(text=element_text(size=14),plot.title = element_text(hjust = 0.5)) +
      labs(x="",y="streamflow (cms)", title=paste0(gage,", ",paste(format(range(dates),"%Y-%m-%d"),collapse=" to ")))

  n_cpd <- max(eventsAll$ix_cpd,na.rm=T)
  if (!is.infinite(n_cpd)) {
    tmp <- subset(eventsAll, !is.na(ix_cpd))
    peaks_cpd <- subset(data2, time %in% tmp$peak_cpd)
    starts_cpd <- subset(data2, time %in% tmp$start_cpd)
    ends_cpd <- subset(data2, time %in% tmp$end_cpd)
    peaks_cpd <- peaks_cpd[!duplicated(peaks_cpd),]
    starts_cpd <- starts_cpd[!duplicated(starts_cpd),]
    ends_cpd <- ends_cpd[!duplicated(ends_cpd),]
    gg1 <- gg1 +
        geom_point(data=peaks_cpd,aes(time,value),color="red",shape=1,size=3) +
        geom_point(data=starts_cpd,aes(time,value),color="blue",shape=1,size=3) +
        geom_point(data=ends_cpd,aes(time,value),color="purple",shape=1,size=3)
  }

  f1 <- paste0("figs/",gage,"_obs_",nwinSpan1,"/event_peak_period_",kk,".png")
  if (!dir.exists(dirname(f1))) dir.create(dirname(f1),recursive=TRUE)
  ggsave(filename=f1,plot=gg1,units="in",width=12.5,height=3.5,dpi=300)

}}}
}
