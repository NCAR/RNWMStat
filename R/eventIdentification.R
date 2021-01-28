# Event identification for hourly streamflow (model or observation)

eventIdentification <- function(data, snowy=FALSE, slow=FALSE,threshPeak,threshold_prob=FALSE,nhourCompound=-1) {

# input arguments
# data: data.frame for streamflow, 1st column is time in POSIXct format,
#       2nd column is the flow values
# snowy: logical, whether the basin is dominated by SEASONAL snow
# slow: logical, whether the basin is a slow-response (e.g.,groundwater-deriven) basin
# threshPeak: threshold for event peaks; peaks below threshold are
#       disgarded
# threshold_prob: logical varible to indicate whether threhPeak
#       is climatological probabilities or the actual values in units of
#       streamflow data
# nhourCompound: max distance (hours) for construct compound events
#       set to a number >=0 if compound events are desired. For example, 
#       if nhourCompound=2, those events next to each other with a distance
#       <=2 hours are combined into a compound event.

# parameters 
nwinSpan=1.5*24
nwinDecay=30*24
minEventDist=24
minRiseDuration=2
minRecessionDuration=6
if (snowy) {
  nwinSpan=10*24
  nwinDecay=90*24
  minRiseDuration=24*5
  minRecessionDuration=24*5
}
if (slow) {
  nwinSpan=10*24
  nwinDecay=60*24
  minRiseDuration=24*2
  minRecessionDuration=24*5
}
minEventDuration=minRiseDuration + minRecessionDuration

data <- na.omit(data)
names(data) <- c("time","value")
if (threshold_prob) {
thresh1 <- quantile(data$value, threshPeak)
} else {
thresh1 <- threshPeak
}
thresh2 <- quantile(data$value, 0.5) #median flow

# identify data gaps and break into chunks with no missing data
dates <- seq(min(data$time),max(data$time), by="hour")
dates1 <- dates[!dates %in% data$time]
chunks <- match(dates1, dates)
nchunk <- length(chunks)+1

# loop through chunks to identy peaks for each chunk and then put them back together
dataAll <- eventsAll <- data.table::data.table()
for (i1 in 1:nchunk) {

   # start index of current non-missing period
   if (i1==1) { j1 <- 1
   } else { j1 <- chunks[i1-1]+1 }

   # end index of current non-missing period
   if (i1==nchunk) { j2 <- length(dates)
   } else { j2 <- chunks[i1]-1 }

   if (j1>j2) next
   if (length(j1:j2) < 6) next 

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
   for (i2 in 3:(nrow(data2)-1))
     if (d1[i2]>=0 & d1[i2+1]<=0) ipeak <- c(ipeak,i2)
   if (length(ipeak)==0) next

   # construct initial events
   events0 <- data.table::data.table()
   for (i2 in 1:length(ipeak)) {

     t2 <- ipeak[i2]

     # end point (on smoothed data)
     ix2 <- t2+nwinDecay
     if (ix2>nrow(data2)) ix2 <- nrow(data2)
     x1 <- data2$smooth[t2:ix2]
     x2 <- diff(x1)
     x3 <- zoo::rollsum(x2,6,align="left")
     ix3 <- which(x3[-(1:minRecessionDuration)] > -0.01)
     if (length(ix3)==0) next
     end1 <- data2$time[t2+min(ix3)+minRecessionDuration]
  
     # start point (on smoothed data)
     ix2 <- t2-nwinDecay
     if (ix2<1) ix2 <- 1
     x1 <- data2$smooth[seq(t2,ix2,-1)]
     x2 <- diff(x1)
     x3 <- zoo::rollsum(x2,6,align="right")
     ix3 <- which(x3[-(1:minRiseDuration)] > -0.01)
     if (length(ix3)==0) next
     start1 <- data2$time[t2-(min(ix3)+minRiseDuration)]
     if (nrow(events0)>=1) 
       if (start1 < events0$end[nrow(events0)]) 
         start1 <- events0$end[nrow(events0)]
  
     # peak point (on original data) 
     dt1 <- subset(data2, time>=start1 & time<=end1)
     peak1 <- dt1$time[which.max(dt1$value)]

     # add to event data frame 
     events0 <- rbind(events0,data.table::data.table(start=start1,peak=peak1,end=end1))
   }
   if (nrow(events0)==0) next

   # compute event duration 
   events0$nhour <- as.integer(difftime(events0$end,events0$start,units="hour"))+1
   events0$nrise <- as.integer(difftime(events0$peak,events0$start,units="hour"))+1
   events0$nrece <- as.integer(difftime(events0$end,events0$peak,units="hour"))

   # remove spurious events
   events0 <- subset(events0, nhour>=6)
   if (nrow(events0)==0) next

   # adjust end points based on start point of next event
   if (nrow(events0)>1) {
   for (i2 in 1:(nrow(events0)-1)) {
     iend <- match(events0$end[i2], data2$time)
     istart <- match(events0$start[i2+1], data2$time)
     if (iend>istart | ((istart-iend)<=10 & data2$value[istart]<data2$value[iend]))
       events0$end[i2] <- data2$time[istart]
   }}

   # discard events at the start or end of periods with data missing
   # if the flow at start/end point is below median flow (thresh2)
   ix1 <- match(events0$start,data2$time)
   ix2 <- which(ix1 != 1 & !is.na(data2$value[ix1-1]))
   ix2 <- ix2 | (ix1==1 & data2$value[ix1]<=thresh2)
   events0 <- events0[ix2,]
   ix1 <- match(events0$end,data2$time)
   ix2 <- which(ix1 != nrow(data2) & !is.na(data2$value[ix1+1]))
   ix2 <- ix2 | (ix1==nrow(data2) & data2$value[ix1]<=thresh2)
   events0 <- events0[ix2,]

   # merge if duplicated start/peak/end
   events0 <- events0[!duplicated(events0),]
   for (tag1 in c("start","peak","end")) {
     dt1 <- subset(data.table::as.data.table(table(events0[[tag1]])),N>1)
     dt1$V1 <- as.POSIXct(dt1$V1, format="%Y-%m-%d %H:%M:%S")
     if (nrow(dt1)>0) {
       for (t1 in dt1$V1) {
         events1 <- subset(events0, get(tag1) == t1)
         events0 <- subset(events0,! get(tag1) %in% t1)
         tmp <- data.table::data.table(start=min(events1$start),
           peak=events1$peak[which.max(data2$value[match(events1$peak,data2$time)])],
           end=max(events1$end))
         tmp$nhour <- as.integer(difftime(tmp$end,tmp$start,units="hour"))+1
         tmp$nrise <- as.integer(difftime(tmp$peak,tmp$start,units="hour"))+1
         tmp$nrece <- as.integer(difftime(tmp$end,tmp$peak,units="hour"))

         events0 <- rbind(events0,tmp) 
   }}}

   # combine events under the following conditions
   # 1. duration too short (rising limb + recession limb),
   #    when nearby event exists (otherwise remove event)
   # 2. peak too close to neighbor event peak
   # 3. event peak falls inside a neighbor event period
   # 4. event too small 
   # 5. one-legged event

   events1 <- data.frame()
   while(nrow(events0)>1) {

   events0 <- events0[order(events0$start),]

   # 1. duration too short
   ix1 <- which(events0$nhour < minEventDuration)
   # 2. peak too close
   ix1 <- c(ix1, which(as.numeric(diff(events0$peak),"hours") < minEventDist)+1)
   # 3. peak inside neighbor event
   ne1 <- nrow(events0)
   peak1 <- events0$peak[1:(ne1-1)]
   start1 <- events0$start[2:ne1]
   peak2 <- events0$peak[2:ne1]
   end2 <- events0$end[1:(ne1-1)]
   ix1 <- c(ix1, which(peak1>=start1), which(peak2<=end2))
   # 4. event too small (based on smoothed data) 
   s2 <- data2$smooth[match(events0$start,data2$time)]
   p2 <- data2$smooth[match(events0$peak,data2$time)]
   e2 <- data2$smooth[match(events0$end,data2$time)]
   h1s <- p2-s2 #rising limb height
   h2s <- p2-e2 #recession limb height
   ix1 <- c(ix1, which(ifelse(h1s<=h2s,h1s,h2s)<0.1))
   # 5. one-legged event (based on original data)
   s1 <- data2$value[match(events0$start,data2$time)]
   p1 <- data2$value[match(events0$peak,data2$time)]
   e1 <- data2$value[match(events0$end,data2$time)]
   h1 <- p1-s1 #rising limb height
   h2 <- p1-e1 #recession limb height
   ix1 <- which(h1<0 | h2<0 | abs(h1)<(0.2*abs(h2)) | abs(h2)<(abs(h1)*0.2))

   if (length(ix1)==0) break
   
   ix1 <- sort(unique(ix1))
   i2 <- ix1[1]
   dif1 <- as.integer(difftime(events0$start[i2],events0$end[i2-1],unit="hours"))
   dif2 <- as.integer(difftime(events0$start[i2+1],events0$end[i2],unit="hours"))
   dist1 <- as.integer(difftime(events0$peak[i2],events0$peak[i2-1],unit="hours"))
   dist2 <- as.integer(difftime(events0$peak[i2+1],events0$peak[i2],unit="hours"))

   flag1 <- 0   
   #first event or event on rising limb, merge with the next event
   if (i2==1 | h2s<0.1 | h2<0 | abs(h2)<(abs(h1)*0.2)) { 
     if(!is.na(dif2) & dif2 < minEventDist) {
       events0$start[i2+1] <- events0$start[i2]; flag1 <- 1
     }
   #last event or event on recession limb, merge with the previous event
   } else if (i2==nrow(events0) | h1s<0.1 | h1<0 | abs(h1)<(abs(h2)*0.5)) {
     if (!is.na(dif1) & dif1 < minEventDist) {
       events0$end[i2-1] <- events0$end[i2]; flag1 <- 1
     } 
   } else {
     if (dist1<=dist2) {
     if (!is.na(dif1) & dif1 < minEventDist) {
       events0$end[i2-1] <- events0$end[i2]; flag1 <- 1
     }} else if (!is.na(dif2) & dif2 < minEventDist) { 
       events0$start[i2+1] <- events0$start[i2]; flag1 <- 1
     }
   }
   if (flag1==0) events1 <- rbind(events1,events0[i2,])
   events0 <- events0[-i2,] # now remove the event
   
   # recompute event duration
   events0$nhour <- as.integer(difftime(events0$end,events0$start,units="hour"))+1
   events0$nrise <- as.integer(difftime(events0$peak,events0$start,units="hour"))+1
   events0$nrece <- as.integer(difftime(events0$end,events0$peak,units="hour"))
   }
  
   # add back those events that are not combined above 
   if (nrow(events1)>0) {
     events0 <- rbind(events0,events1)
     events0 <- events0[order(events0$start),]
   } 

   if (nrow(events0)>1) {
   # adjust peaks, start and end points
   # start point should be the lowerest point on the rising limb
   # end point should be the lowerest point on the recession limb
   # peak should be the highest point during event period
   for (i2 in 1:nrow(events0)) {
     istart <- match(events0$start[i2], data2$time)
     ipeak <- match(events0$peak[i2], data2$time)
     iend <- match(events0$end[i2], data2$time)
   
     ix0 <- which.min(data2$value[istart:ipeak])
     events0$start[i2] <- data2$time[istart+ix0-1]
   
     ix0 <- which.min(data2$value[ipeak:iend])
     events0$end[i2] <- data2$time[ipeak+ix0-1]
   
     dt1 <- subset(data2, time %in% seq(events0$start[i2],events0$end[i2],by="hour"))
     events0$peak[i2] <- dt1$time[which.max(dt1$value)]
   }
  
   # check if event start/peak/end are in order
   events0 <- events0[order(events0$start),]
   ne1 <- nrow(events0)
   ix1 <- which((1:ne1) != order(events0$peak))
   if (length(ix1)>0) print("WARNING: events peak not in order!")
   ix1 <- which((1:ne1) != order(events0$end))
   if (length(ix1)>0) print("WARNING: events end not in order!")
   
   ix1 <- as.integer(difftime(events0$start[2:ne1],events0$end[1:(ne1-1)],units="hour"))
   if(sum(ix1<0)) print(paste0("WARNING: event starts before previous event ends: ",paste(which(ix1<0)+1,collapse=", ")))

   # combine adjacent events (to form compound events)
   if (nhourCompound >= 0) {
   ne1 <- nrow(events0)
   ends1 <- events0$end[1:(ne1-1)]
   starts1 <- events0$start[2:ne1]
   dif1 <- as.integer(difftime(starts1,ends1,units="hour"))
   ix1 <- which(dif1<=nhourCompound)+1
   if (length(ix1)==0) break
   kk <- 1
   while(kk<=length(ix1)) {
     ixs <- ix1[kk]-1
     while(1) {
       ixs <- c(ixs,ix1[kk])
       kk <- kk +1
       if (kk>length(ix1)) break
       if (ix1[kk] > (ix1[kk-1]+1)) break
     }
     peaks <- subset(data2, time %in% events0$peak[ixs])$value
     events0$peak[ixs[1]] <- events0$peak[ixs][which.max(peaks)]
     events0$end[ixs[1]] <- events0$end[ixs[length(ixs)]]
   }
   events0 <- events0[!(1:ne1 %in% ix1),]
   }

   # adjust adjecent events so that the previous event ends at lowest point
   # between the two peaks and the next event starts at the same point
   for (i2 in 1:(nrow(events0)-1)) {
     ipeak <- match(events0$peak[i2],data2$time)
     ipeak1 <- match(events0$peak[i2+1],data2$time)
     values <- data2$value[ipeak:ipeak1]
     if (sum(values>=thresh2)==length(values)) {
       ix0 <- which.min(values)
       events0$end[i2] <- events0$start[i2+1] <- data2$time[ipeak+ix0-1]
   }}

   }

   events0$nhour <- as.integer(difftime(events0$end,events0$start,units="hour"))+1
   events0$nrise <- as.integer(difftime(events0$peak,events0$start,units="hour"))+1
   events0$nrece <- as.integer(difftime(events0$end,events0$peak,units="hour"))

   # remove events below threshold
   ix1 <- which(data2$value[match(events0$peak,data2$time)] >= thresh1)
   events0 <- events0[ix1,]

   # combine events and data from all chunks
   data2$hour <- NULL
   dataAll <- rbind(dataAll, data2)
   eventsAll <- rbind(eventsAll, events0)
}

# put back those chunks that have no peaks (and hence removed during event identification)
data1 <- subset(data, ! time %in% dataAll$time)
if (nrow(data1)>0) {
data1$smooth <- NA
dataAll <- rbind(dataAll, data1)
dataAll <- dataAll[order(dataAll$time),]
}

return(list(eventsAll=eventsAll, dataAll=dataAll))
}

