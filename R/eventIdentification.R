eventIdentification <- function(data,tag,nwinSpan=36,plot=TRUE,
                                probPeak=0.9,probLowflow=0.5,probRange=0.3,
                                minEventDuration=6,maxRiseDuration=120,maxRecessionDuration=240,
                                nwinShift=12,minInterval=6,minLengthData=6) {

  ########### definition of parameters #####################
  # data: data frame for streamflow, 1st column is time in POSIXct format,
  #       2nd column is the flow values
  # tag: string to name the data set, e.g., "14400000_mod","14400000_obs"
  # nwinSpan: window size (in hours) for smoothing
  # nwinShift: size of window (nwinShift*2+1) around smoothed peak to identify actual peak,
  #       as there often exists a shift between smooth and actual peaks (hours)
  # minInterval: minimum seperation between peaks; peaks that are too close
  #       will be combined
  # minEventDuration: minimum event duration in hours; events with duration
  #       shorter than the threshold will be combined with neighbour events
  # maxRiseDuration: maximum duration in hours for the rise limb
  # maxRecessionDuration: maximum duration in hours for the recession limb
  # probPeak: probability threshold for event peaks; peaks below threshold are
  #       disgarded
  # probLowflow: probability threshold for low flows, used to determine event
  #       start/end points
  # probFlowRange: probability threshold for event flow range (i.e., flow
  #       difference between peak and start point, or between peak and end point).
  #       Events with flow ranges below this threshold are discarded
  # minLengthData: mininum length of record to perform event separation (hours)
  #       flow time series is first devided into chunks with no missing values;
  #       chunks too short are discarded

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
    iStarts0 <- sapply(iPeaks, function(x) max(iLowflows[iLowflows<x]))
    iEnds0 <- sapply(iPeaks, function(x) min(iLowflows[iLowflows>x]))
    iStarts0[is.infinite(iStarts0) | is.na(iEnds0)] <- 1
    iEnds0[is.infinite(iEnds0) | is.na(iEnds0)] <- nrow(df1)

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
    eventDt <- data.table::data.table(start=df1$time[iStarts],
                          peak=df1$time[iPeaks],
                          end=df1$time[iEnds])

    # fine tune start and end points
    eventDt <- tuneEventStartEnd(eventDt, df1)

    eventDt
  }

  ############## processing starts here #############
  data <- na.omit(data)
  names(data) <- c("time","value")
  thresh1 <- quantile(data$value, probPeak)
  thresh2 <- quantile(data$value, probLowflow)
  thresh3 <- quantile(data$value, probRange)

  # minor adjustments to peak threshold
  if (thresh1 <= max(data$value)*0.01) thresh1 <- max(data$value)*0.01
  if (thresh1 <= 1.0) thresh1 <- 1.0 #cms, convert to cfs if flow unit is cfs

  # identify data gaps and break into chunks with no missing data
  dates <- seq(min(data$time),max(data$time), by="hour")
  dates1 <- dates[!dates %in% data$time]
  chunks <- match(dates1, dates)
  nchunk <- length(chunks)+1

  # loop through chunks to identy peaks for each chunk and then put them back together
  dataAll <- data.table::data.table()
  eventsAll <- data.table::data.table()
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

    data2$hour <- 1:nrow(data2)

    #kernel smoothing
    #fit <- with(data2,ksmooth(hour, value, kernel="normal",bandwidth=span))
    #data2$smooth <- fit$y

    #local weighted regression
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
    peaks1 <- data.table::data.table()
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

    # remove peaks with small peak-valley difference
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

    #print(paste0(i1," ", nrow(events1)," ",nrow(events2)," ",nrow(events3)))

    # put them back together
    data2$hour <- NULL
    dataAll <- rbind(dataAll, data2)
    eventsAll <- rbind(eventsAll, events3)
  }


  outfile <- paste0("events/events_",tag,".Rdata")
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile),recursive=TRUE)
  save(eventsAll,file=outfile)

  if (plot) {

    # add in those chunks that have no peaks
    data1 <- subset(data, ! time %in% dataAll$time)
    data1$smooth <- NA
    dataAll <- rbind(dataAll, data1)
    dataAll <- dataAll[order(time),]

    #plot hydrographs and peaks (in half-year chunks)
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
        gg1 <- ggplot2::ggplot(data=data2,aes(time,value)) +
          ggplot2::geom_line(color="darkgrey") +
          ggplot2::geom_line(aes(time,smooth),color="black",size=0.4)+
          ggplot2::geom_point(data=peaks,aes(time,value),color="red",shape=2,size=1) +
          ggplot2::geom_point(data=starts,aes(time,value),color="blue",shape=8,size=1) +
          ggplot2::geom_point(data=ends,aes(time,value),color="green",shape=1,size=1) +
          ggplot2::geom_hline(yintercept=thresh1,linetype="dashed",color="black",size=0.5) +
          ggplot2::theme(text=element_text(size=14),plot.title = element_text(hjust = 0.5)) +
          ggplot2::labs(x="",y="streamflow (cms)", title=paste0(gage,", ",paste(format(range(dates),"%Y-%m-%d"),collapse=" to ")))

        f1 <- paste0("figs/",tag,"_",nwinSpan,"/event_peak_",tag,"_",kk,".png")
        if (!dir.exists(dirname(f1))) dir.create(dirname(f1),recursive=TRUE)
        ggplot2::ggsave(filename=f1,plot=gg1,units="in",width=12.5,height=3.5,dpi=300)

      }}}

  eventsAll
}
