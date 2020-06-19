# matched observed events with modelled events

matchEvents <- function(eventsAllMod, eventsCompoundMod,
  eventsCompoundObs,maxDist=48) {

###### Inputs ##############
# eventsAllMod: data frame listing all single-peak events for model
#    streamflow (identified by eventIdentification.R)
# eventsCompoundMod: data frame of compound events for model streamflow
# eventsCompoundObs: data frame of compound events for observed streamflow
# maxDist: window size (in hours) for identifying model events for
#    a given observed event. Consider starting with 48 hours for
#    non-snow basins and 7*24 (i.e., one week) for snow basins

##### Output ##############
# data frame with infomation on events matching. It has the same
# number of rows as eventsCompoundObs
#   1st column: index of observed events (in eventsCompoundObs)
#   2nd column: match category 
#       1 = matched with eventsCompoundMod
#       2 = matched with eventsAllMod
#       3 = missed by model
#   3rd column: index of matched model events in eventsCompoundMod
#   4th column: index of matched model events in eventsAllMod
 

no1 <- nrow(eventsCompoundObs)
nm1 <- nrow(eventsCompoundMod)
nm2 <- nrow(eventsAllMod) 

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

dfMatch

}
