# Event identification for hourly streamflow (model or observation)
 
eventIdentification <- function(data,nwinSpan=36,threshold_prob=TRUE,
  threshPeak=0.9,threshLowflow=0.5,threshFlowRange=0.3,maxGapFill=5,
  minEventDuration=6,maxRiseDuration=120,maxRecessionDuration=240,
  maxPeakDistCompound=120, nwinShift=12,minInterval=6,minLengthData=6) {

########### Inputs #####################
# data: data.frame for streamflow, 1st column is time in POSIXct format, 
#       2nd column is the flow values
# nwinSpan: window size (in hours) for smoothing; consider starting with
#       240 hours for snow-driven streamflow with diel fluctuations and
#       36 hours for non-snow-driven streamflow  
# nwinShift: size of window (nwinShift*2+1) around smoothed peak to 
#       identify actual peak, as there often exists a shift between smoothed
#       and actual peaks (hours)
# minInterval: minimum seperation between peaks; peaks that are too close 
#       will be combined
# minEventDuration: minimum event duration in hours; events with duration 
#       shorter than the threshold will be combined with neighbour events
# maxRiseDuration: maximum duration in hours for the rise limb
# maxGapFill: maximum duration (in hours) of data gaps to be filled with 
#       spline interpolation; gaps larger than this threshold will not be filled
# maxRecessionDuration: maximum duration in hours for the recession limb
# maxPeakDistCompound: maximum distance between adjacent peaks that can be
#       considered combined into compound events
# threshold_prob: logical varible to indicate whether the three threhold 
#       parametes (threshPeak, threshLowflow, threshFlowRange) are probablistic
#       or not
# threshPeak: threshold for event peaks; peaks below threshold are 
#       disgarded
# threshLowflow: threshold for low flows, used to determine event 
#       start/end points
# threshFlowRange: probability threshold for event flow range (i.e., flow 
#       difference between peak and start point, or between peak and end point).
#       Events with flow ranges below this threshold are discarded
# minLengthData: mininum length of record to perform event separation (hours)
#       flow time series is first devided into chunks with no missing values;
#       chunks too short are discarded 

######## outputs ##############
# Output is a list of three data frames
# 1st data frame (eventsAll): all single-peak events identified, with three
#       POSIXct columns indicating the start, peak, and end hours of all events
# 2nd data frame (eventsCompound): final list of events identified where adjcent 
#       single-peak events are combined into a compound event
# 3rd data frame: the original streamflow data frame with short gap filled 
#       with spline interpolation, as well as the smoothed time series   


# function to get the index of the last minimum value
# which.min gets the first minimum value index
which.min2 <- function(x, last.index = FALSE, ...){
  if(last.index) max(which(x == min(x, ...))) else which.min(x)
}

# function to fine tune start and end points
# start point should have the minimum flow value on the rising limb
# end point should have the minimum flow value on the recession limb
tuneEventStartEnd <- function(events, df1) {
   for (k1 in 1:nrow(events)) {

     iStart <- match(events$start[k1],df1$time)
     iPeak <- match(events$peak[k1],df1$time)
     iEnd <- match(events$end[k1],df1$time)

     idx1 <- which.min2(df1$value[iStart:iPeak],last.index=TRUE)
     if (idx1!=1) events$start[k1] <- events$start[k1]+(idx1-1)*3600

     idx1 <- which.min(df1$value[iPeak:iEnd])
     if (idx1!=(iEnd-iPeak+1)) events$end[k1] <- events$peak[k1]+(idx1-1)*3600
   }
   events
}

# function to identify event start/end points given peak time and data
getEvents <- function(tPeaks,df1,thresh) {

iLowflows <- which(df1$value < thresh)
iPeaks <- match(tPeaks,df1$time)
iStarts0 <- sapply(iPeaks, function(x) max(iLowflows[iLowflows<x]-1))
iEnds0 <- sapply(iPeaks, function(x) min(iLowflows[iLowflows>x]+1))
iStarts0[is.infinite(iStarts0) | is.na(iEnds0) | iStarts0<1] <- 1
iEnds0[is.infinite(iEnds0) | is.na(iEnds0) | iEnds0>nrow(df1)] <- nrow(df1)

n1 <- length(iPeaks)
iStarts1 <- iEnds1 <- rep(NA,n1)
for (k1 in 1:n1) {
  if (k1 == 1) {
    iStarts1[k1] <- iStarts0[k1]
    if (n1==1) iEnds1[k1] <- iEnds0[k1]
  } else {
    j1 <- iPeaks[k1-1]-1+which.min(df1$value[iPeaks[k1-1]:iPeaks[k1]])
    iStarts1[k1] <- j1
    iEnds1[k1-1] <- j1
    if (k1 == n1) iEnds1[k1] <- iEnds0[k1]
}}

iStarts <- sapply(1:n1, function(x) max(iStarts0[x],iStarts1[x]))
iEnds <- sapply(1:n1, function(x) min(iEnds0[x],iEnds1[x]))
iPeaks <- sapply(1:length(iPeaks),function(x) which.max(df1$value[iStarts[x]:iEnds[x]])+iStarts[x]-1)

# events (start time, peak time, end time)
eventDt <- data.frame(start=df1$time[iStarts],
                      peak=df1$time[iPeaks],
                      end=df1$time[iEnds])

# fine tune start and end points 
eventDt <- tuneEventStartEnd(eventDt, df1)

eventDt
}

# function to identify compound events given a series single events
compoundEvents <- function(events,df1,thresh,ix_cpd0) {

ne1 <- nrow(events)
events$ix_cpd <- NA
events$start_cpd <- events$start
events$peak_cpd <- events$peak
events$end_cpd <- events$end

# first label all the single events that belong to compound events
if (ne1>1) {
kk <- ix_cpd0
for (k1 in 1:(ne1-1)) 
  if (as.integer(difftime(events$start[k1+1],events$end[k1],units="hours"))<=5) {
  if (as.integer(difftime(events$peak[k1+1],events$peak[k1],units="hours"))<=maxPeakDistCompound) {
    min0 <- min(subset(df1,time %in% seq(events$end[k1],events$start[k1+1],by="hour"))$value,na.rm=T) 
    if (min0 > thresh) {
      if (!is.na(events$ix_cpd[k1])) {
        events$ix_cpd[k1+1] <- events$ix_cpd[k1]
      } else {
        kk <- kk + 1; events$ix_cpd[k1:(k1+1)] <- kk
      }}}}

# then identify the start, peak, and end points of compound events
if (kk>ix_cpd0) {
  for (k1 in (ix_cpd0+1):kk) {
    
     ix1 <- which(!is.na(events$ix_cpd) & events$ix_cpd==k1)
     events$start_cpd[ix1] <- events$start[min(ix1)]
     events$end_cpd[ix1] <- events$end[max(ix1)]
     ipeaks <- events$peak[ix1]
     events$peak_cpd[ix1] <- ipeaks[which.max(df1$value[match(ipeaks,df1$time)])]

}}
}
events
}

############## processing starts here #############
data <- na.omit(data)
names(data) <- c("time","value")
if (threshold_prob) {
thresh1 <- quantile(data$value, threshPeak)
thresh2 <- quantile(data$value, threshLowflow)
thresh3 <- quantile(data$value, threshFlowRange)
} else {
thresh1 <- threshPeak
thresh2 <- threshLowflow
thresh3 <- threshFlowRange
}

# minor adjustments to peak threshold
if (thresh1 <= max(data$value)*0.01) thresh1 <- max(data$value)*0.01
if (thresh1 <= 1.0) thresh1 <- 1.0 #cms, convert to cfs if flow unit is cfs

# fill short data gaps with spline interpolation
dates <- seq(min(data$time),max(data$time), by="hour")
dates1 <- dates[!dates %in% data$time]
if (length(dates1)>0) {
data <- rbind(data,data.frame(time=dates1,value=NA))
data <- data[order(data$time),]
data$value <- zoo::na.approx(data$value,maxgap=maxGapFill,na.rm=FALSE)
data <- subset(data,!is.na(value))
}

# identify remaining data gaps and break into chunks with no missing data
dates <- seq(min(data$time),max(data$time), by="hour")
dates1 <- dates[!dates %in% data$time]
chunks <- match(dates1, dates)
nchunk <- length(chunks)+1

# loop through chunks to identy peaks for each chunk and then put them back together
dataAll <- eventsAll <- data.frame()
for (i1 in 1:nchunk) {

   # start index of current non-missing period
   if (i1==1) { j1 <- 1
   } else { j1 <- chunks[i1-1]+1 }

   # end index of current non-missing period
   if (i1==nchunk) { j2 <- length(dates)
   } else { j2 <- chunks[i1]-1 }

   if (j1>j2) next
   if (length(j1:j2) < minLengthData) next 

   # data for current chunk
   data2 <- subset(data, time %in% dates[j1:j2])
   if (max(data2$value) < thresh1) next

   #local weighted regression smoothing
   data2$hour <- 1:nrow(data2)
   span <- nwinSpan/nrow(data2)
   fit <- loess(value ~ hour, degree=1,span = span, data=data2)
   data2$smooth <- fit$fitted

   # identify peaks in smoothed data
   d1 <- c(NA, diff(data2$smooth))
   ipeak <- NULL
   for (i2 in 2:(nrow(data2)-1)) 
     if (d1[i2]>=0 & d1[i2+1]<=0) ipeak <- c(ipeak,i2)
   if (length(ipeak)==0) next

   # identify corresponding peaks in the original data
   ipeak1 <- rep(NA, length(ipeak))
   for (i2 in 1:length(ipeak)) {
     j1 <- ipeak[i2]-nwinShift
     j2 <- ipeak[i2]+nwinShift
     if (j1<1) j1 <- 1
     if (j2>nrow(data2)) j2 <- nrow(data2)
     ix2 <- which.max(data2$value[j1:j2])
     ipeak1[i2] <- ipeak[i2]-nwinShift-1+ix2
   }
   ipeak1 <- unique(ipeak1)
   ipeak1[ipeak1 < 1] <- 1
   ipeak1[ipeak1 > nrow(data2)] <- nrow(data2)

   # peaks identified
   peaks <- data2[ipeak1,]

   # remove those below the threshold
   peaks <- subset(peaks, value >= thresh1)
   if (nrow(peaks)==0) next

   # combine those peaks that are too close
   int1 <- which(as.numeric(diff(peaks$time),"hours") < minInterval)
   rowIdx <- 1:nrow(peaks)
   rowIdx <- rowIdx[!rowIdx %in% (int1+1)]
   peaks1 <- data.frame()
   for (k1 in rowIdx) {
      if (k1 %in% int1) {
        if (peaks$value[k1] > peaks$value[k1+1]) {
          peaks1 <- rbind(peaks1,peaks[k1,])
        } else {
          peaks1 <- rbind(peaks1,peaks[k1+1,])
        }
      } else {
        peaks1 <- rbind(peaks1,peaks[k1,])
      }
   }
   if (nrow(peaks1)==0) next

   # identify start and end points of all events
   events1 <- getEvents(peaks1$time, data2,thresh2)

   # remove events with a rising/recession limb that is too short (verticallly)
   iPeaks <- match(events1$peak,data2$time)
   iStarts <- match(events1$start,data2$time)
   iEnds <- match(events1$end,data2$time)
   ix1 <- NULL
   n1 <- length(iPeaks)
   for (k1 in 1:n1)
     if (min(data2$value[iPeaks[k1]]-data2$value[iStarts[k1]],data2$value[iPeaks[k1]]-data2$value[iEnds[k1]]) >= thresh3) ix1 <- c(ix1,k1)
   if (length(ix1)==0) next
   events2 <- getEvents(events1$peak[ix1], data2,thresh2)

   # remove events that are too short in duration
   ix1 <- which(as.integer(difftime(events2$end,events2$start,units="hours")) >= minEventDuration)
   if (length(ix1)==0) next
   events3 <- getEvents(events2$peak[ix1], data2,thresh2)

   # adjust long starts and long tails
   for (k1 in 1:nrow(events3)) {
      t1 <- as.integer(difftime(events3$peak[k1],events3$start[k1],units="hours"))
      if (t1>maxRiseDuration) events3$start[k1] <- events3$peak[k1]-maxRiseDuration*3600
      t1 <- as.integer(difftime(events3$end[k1],events3$peak[k1],units="hours"))
      if (t1>maxRecessionDuration) events3$end[k1] <- events3$peak[k1]+maxRecessionDuration*3600
  }

  # final adjustments to start and end points
  events3 <- tuneEventStartEnd(events3,data2)

  # identify compound events
  ix_cpd0 <- 0
  if (nrow(eventsAll)>=1) ix_cpd0 <- max(eventsAll$ix_cpd,na.rm=T)
  if (is.infinite(ix_cpd0)) ix_cpd0 <- 0
  events3 <- compoundEvents(events3, data2, thresh2,ix_cpd0)

  # put them back together
  data2$hour <- NULL
  dataAll <- rbind(dataAll, data2)
  eventsAll <- rbind(eventsAll, events3)
}

# add in those chunks that have no peaks
data1 <- subset(data, ! time %in% dataAll$time)
if (nrow(data1)>0) {
data1$smooth <- NA
dataAll <- rbind(dataAll, data1)
dataAll <- dataAll[order(dataAll$time),]
}

# compile event list that combines regular events with compound events
eventsAll$start0 <- eventsAll$peak0 <- eventsAll$end0 <- eventsAll$end
ix1 <- which(is.na(eventsAll$ix_cpd))
eventsAll$start0[ix1] <- eventsAll$start[ix1]
eventsAll$peak0[ix1] <- eventsAll$peak[ix1]
eventsAll$end0[ix1] <- eventsAll$end[ix1]
ix1 <- which(!is.na(eventsAll$ix_cpd))
eventsAll$start0[ix1] <- eventsAll$start_cpd[ix1]
eventsAll$peak0[ix1] <- eventsAll$peak_cpd[ix1]
eventsAll$end0[ix1] <- eventsAll$end_cpd[ix1]

eventsCompound <- eventsAll[,c("start0","peak0","end0")]
eventsCompound <- eventsCompound[!duplicated(eventsCompound),]
names(eventsCompound) <- c("start","peak","end")

eventsAll <- eventsAll[,c("start","peak","end")]

# discard events at the start or end of periods with data missing
ix1 <- match(eventsAll$start,dataAll$time)
ix2 <- which(ix1 != 1 & !is.na(dataAll$value[ix1-1]))
ix2 <- c(ix2,which(ix1==1 & dataAll$value[ix1]<=thresh2))
eventsAll <- eventsAll[ix2,]
ix1 <- match(eventsAll$end,dataAll$time)
ix2 <- c(ix2, which(ix1==nrow(dataAll) & dataAll$value[ix1]<=thresh2))
ix2 <- which(ix1 != nrow(dataAll) & !is.na(dataAll$value[ix1+1]))
eventsAll <- eventsAll[ix2,]

ix1 <- match(eventsCompound$start,dataAll$time)
ix2 <- which(ix1 != 1 & !is.na(dataAll$value[ix1-1]))
ix2 <- c(ix2,which(ix1==1 & dataAll$value[ix1]<=thresh2))
eventsCompound <- eventsCompound[ix2,]
ix1 <- match(eventsCompound$end,dataAll$time)
ix2 <- which(ix1 != nrow(dataAll) & !is.na(dataAll$value[ix1+1]))
ix2 <- c(ix2, which(ix1==nrow(dataAll) & dataAll$value[ix1]<=thresh2))
eventsCompound <- eventsCompound[ix2,]

list(eventsAll,eventsCompound,dataAll)

}
