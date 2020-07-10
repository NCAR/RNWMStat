# Match observed (compound) events with model events

matchEvents <- function(data_mod, eventsAllMod, eventsCompoundMod,
  eventsCompoundObs,maxDist=48) {

###### Inputs ##############
# data_mod: data.frame for model streamflow, 1st column is time in POSIXct format,
#       2nd column is the flow values
# eventsAllMod: data frame listing all single-peak events for model
#    streamflow (identified by eventIdentification.R)
# eventsCompoundMod: data frame of compound events for model streamflow
# eventsCompoundObs: data frame of compound events for observed streamflow
# maxDist: window size (in hours) for identifying model events for
#    a given observed event. Consider starting with 48 hours for
#    non-snow basins and 7*24 (i.e., one week) for snow basins

##### Output ##############
# data frame with columns corresponding to the start, peak, and end times of
# matched observed and model events. It has the same number of rows as 
#   eventsCompoundObs
#   1st - 3rd columns: start, peak and end times of observed events
#   4th - 6th columns: start, peak and end times of model events 
#   7th column: match category 
#       1 = matched with eventsCompoundMod
#       2 = matched with eventsAllMod
#       3 = missed by model

no1 <- nrow(eventsCompoundObs)
nm1 <- nrow(eventsCompoundMod)
nm2 <- nrow(eventsAllMod) 

if (no1==0 | (nm1+nm2)==0) { 
  print("WARNING: no events to match")
  return(data.frame())
}
 
dfMatch <- data.frame(ix_obs=1:no1,match=rep(NA,no1),
  ix_mod1=rep(NA,no1), ix_mod2=rep(NA,no1))

for (k1 in 1:no1) {

  # first round: match observed compound events with model compound events
  idx1 <- which(! (1:nm1) %in% dfMatch$ix_mod1[k1])
  dist1 <- abs(as.integer(difftime(eventsCompoundObs$peak[k1],eventsCompoundMod$peak[idx1],units="hours")))
  mdist1 <- min(dist1)
  i1 <- idx1[which.min(dist1)]
  if (mdist1 <= maxDist) {
     dfMatch$match[k1] <- 1
     dfMatch$ix_mod1[k1] <- i1
     next
  } else {
  # 2nd round: match remaining observed events with model single-peak events
    idx1 <- which(! (1:nm2) %in% dfMatch$ix_mod2[k1])
    dist1 <- abs(as.integer(difftime(eventsCompoundObs$peak[k1],eventsAllMod$peak[idx1],units="hours")))
    mdist1 <- min(dist1)
    i1 <- idx1[which.min(dist1)]
    if (mdist1 <= maxDist) {
       dfMatch$match[k1] <- 2
       dfMatch$ix_mod2[k1] <- i1
       next
    } else {
   # remaining observed events that are unmatched (i.e., missed by model)
       dfMatch$match[k1] <- 3
    }
  }
}

# assemble matched events of obs and mod
peak_obs <- eventsCompoundObs$peak
start_obs <- eventsCompoundObs$start
end_obs <- eventsCompoundObs$end

peak_mod <- peak_obs
start_mod <- start_obs
end_mod <- end_obs

ix1 <- which(dfMatch$match==1)
if (length(ix1)>0) {
  peak_mod[ix1] <- eventsCompoundMod$peak[dfMatch$ix_mod1[ix1]]
  start_mod[ix1] <- eventsCompoundMod$start[dfMatch$ix_mod1[ix1]]
  end_mod[ix1] <- eventsCompoundMod$end[dfMatch$ix_mod1[ix1]]
}
ix1 <- which(dfMatch$match==2)
if (length(ix1)>0) {
  peak_mod[ix1] <- eventsAllMod$peak[dfMatch$ix_mod2[ix1]]
  start_mod[ix1] <- eventsAllMod$start[dfMatch$ix_mod2[ix1]]
  end_mod[ix1] <- eventsAllMod$end[dfMatch$ix_mod2[ix1]]
}
ix1 <- which(dfMatch$match==3)
if (length(ix1)>0) {
names(data_mod) <- c("time","value")
for (k1 in ix1) {
tmp <- subset(data_mod,time>=eventsCompoundObs$start[k1] & time<=eventsCompoundObs$end[k1])
peak_mod[k1] <- tmp$time[which.max(tmp$value)]
}}

eventsMatched <- data.frame(start_obs,peak_obs,end_obs,start_mod,peak_mod,end_mod)
eventsMatched$match <- dfMatch$match
eventsMatched

}
