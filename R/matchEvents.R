# Match observed events with model events
# Output is a list: 
#     events_obs_match - observed events that are matched with model events
#     events_mod_match - model events that are matched with obs events
#     matched_mod - indication of whether the orginal model event is matched
#     matched_obs - indication of whether the orginal obs event is matched

matchEvents <- function(data_obs,data_mod,events_obs, events_mod, threshold1,snowy=FALSE) {

win1 <- 5*24 # window to look for matching events
if (snowy) win1 <- 3*30*24

no1 <- nrow(events_obs)
nm1 <- nrow(events_mod)
events_obs <- events_obs[order(events_obs$start),]
events_mod <- events_mod[order(events_mod$start),]
events_mod1 <- data.table::data.table()
match_mod <- rep(0,nm1)
match_obs <- rep(0,no1)

# loop through observed events to find matches in the model events
for (i1 in 1:no1) {

dates1 <- seq(events_obs$start[i1],events_obs$end[i1],by="hour")
dates2 <- seq(min(dates1)-win1*3600, max(dates1)+win1*3600,by="hour")
ix0 <- which(events_mod$peak %in% dates2) 
ix1 <- NULL
for (i2 in ix0) {
  dates <- seq(events_mod$start[i2],events_mod$end[i2],by="hour")
  if (sum(dates %in% dates1)>=(min(length(dates1),length(dates))*0.2)) ix1 <- c(ix1,i2)
}

if (length(ix1)==0) { #no matching model events
  d1 <- subset(data_mod, time %in% dates1)
  peak_time <- d1$time[which.max(d1$value)]
  peak_value <- max(d1$value)
  if (peak_value >= threshold1) match_obs[i1] <- 1 
  t1 <- data.table::data.table(start=events_obs$start[i1],
    peak=peak_time, end=events_obs$end[i1])
  t1$nhour <- as.integer(difftime(t1$end,t1$start,units="hour"))+1
  t1$nrise <- as.integer(difftime(t1$peak,t1$start,units="hour"))+1
  t1$nrece <- as.integer(difftime(t1$end,t1$peak,units="hour"))
  events_mod1 <- rbind(events_mod1, t1)
     
} else if (length(ix1)==1) { #find one matching model event
  events_mod1 <- rbind(events_mod1,events_mod[ix1,])
  match_mod[ix1] <- 1
  match_obs[i1] <- 1

} else { 
  # find mulitple model events within observed event period
  # merge them into one compound mode event
  e1 <- events_mod[ix1,]
  peak_time <- e1$peak[which.max(data_mod$value[match(e1$peak,data_mod$time)])]
  t1 <- data.table::data.table(start=min(e1$start),
    peak=peak_time, end=max(e1$end))
  t1$nhour <- as.integer(difftime(t1$end,t1$start,units="hour"))+1
  t1$nrise <- as.integer(difftime(t1$peak,t1$start,units="hour"))+1
  t1$nrece <- as.integer(difftime(t1$end,t1$peak,units="hour"))
  events_mod1 <- rbind(events_mod1, t1)
  match_mod[ix1] <- 1
  match_obs[i1] <- 1
  
}
} # loop events

# make sure events are not overlapping each other
# first model events
n1 <- nrow(events_mod1)
if (n1>1) {
while(1) {
ix0 <- which(duplicated(events_mod1))
events <- events_mod1[ix0,]
events <- events[!duplicated(events),]
events1 <- rbind(events,events_mod1)
jx0 <- which(duplicated(events1))
jx0 <- jx0-nrow(events)

ix1 <- which(events_mod1$end[1:(n1-1)] > events_mod1$start[2:n1])
ix1 <- ix1[!ix1 %in% jx0]
if (length(ix1)==0) break

j1 <- ix1[1]; j2 <- ix1[1]+1
if (events_mod1$start[j2]<=events_mod1$start[j1] & events_mod1$end[j2]>=events_mod1$end[j1]) {
events_mod1[j1,] <- events_mod1[j2,]
} else if (events_mod1$start[j1]<=events_mod1$start[j2] & events_mod1$end[j1]>=events_mod1$end[j2]) {
events_mod1[j2,] <- events_mod1[j1,]
} else {
dates <- seq(events_mod1$start[ix1[1]+1],events_mod1$end[ix1[1]],by="hour")
values <-  data_mod$value[match(dates,data_mod$time)]
events_mod1$start[ix1[1]+1] <- events_mod1$start[ix1[1]+1] + (which.min(values)-1)*3600
events_mod1$end[ix1[1]] <- events_mod1$start[ix1[1]+1]
}}
events_mod1$nhour <- as.integer(difftime(events_mod1$end,events_mod1$start,units="hour"))+1
events_mod1$nrise <- as.integer(difftime(events_mod1$peak,events_mod1$start,units="hour"))+1
events_mod1$nrece <- as.integer(difftime(events_mod1$end,events_mod1$peak,units="hour"))

# then observed events
events_obs1 <- events_obs
while(1) {
ix0 <- which(duplicated(events_obs1))
events <- events_obs1[ix0,]
events <- events[!duplicated(events),]
events1 <- rbind(events,events_obs1)
jx0 <- which(duplicated(events1))
jx0 <- jx0-nrow(events)

ix1 <- which(events_obs1$end[1:(n1-1)] > events_obs1$start[2:n1])
ix1 <- ix1[!ix1 %in% jx0]
if (length(ix1)==0) break
j1 <- ix1[1]; j2 <- ix1[1]+1
if (events_obs1$start[j2]<=events_obs1$start[j1] & events_obs1$end[j2]>=events_obs1$end[j1]) {
events_obs1[j2,] <- events_obs1[j1,]
} else if (events_obs1$start[j1]<=events_obs1$start[j2] & events_obs1$end[j1]>=events_obs1$end[j2]) {
events_obs1[j1,] <- events_obs1[j2,]
} else {
dates <- seq(events_obs1$start[ix1[1]+1],events_obs1$end[ix1[1]],by="hour")
values <-  data_obs$value[match(dates,data_obs$time)]
events_obs1$start[ix1[1]+1] <- events_obs1$start[ix1[1]+1] + (which.min(values)-1)*3600
events_obs1$end[ix1[1]] <- events_obs1$start[ix1[1]+1]
}}
events_obs1$nhour <- as.integer(difftime(events_obs1$end,events_obs1$start,units="hour"))+1
events_obs1$nrise <- as.integer(difftime(events_obs1$peak,events_obs1$start,units="hour"))+1
events_obs1$nrece <- as.integer(difftime(events_obs1$end,events_obs1$peak,units="hour"))


# if duplicated model events, combine obs events
while(1) {
ix1 <- which(duplicated(events_mod1))
if (length(ix1)==0) break

events_obs1$end[ix1[1]-1] <- events_obs1$end[ix1[1]]
peak0 <- data_obs$value[match(events_obs1$peak[ix1[1]-1],data_obs$time)]
peak1 <- data_obs$value[match(events_obs1$peak[ix1[1]],data_obs$time)]
if (peak0<peak1) events_obs1$peak[ix1[1]-1] <- events_obs1$peak[ix1[1]]

events_obs1 <- events_obs1[-ix1[1],]
events_mod1 <- events_mod1[-ix1[1],]
}

# if duplicated obs events, combine model events
while(1) {
ix1 <- which(duplicated(events_obs1))
if (length(ix1)==0) break

events_mod1$end[ix1[1]-1] <- events_mod1$end[ix1[1]]
peak0 <- data_mod$value[match(events_mod1$peak[ix1[1]-1],data_mod$time)]
peak1 <- data_mod$value[match(events_mod1$peak[ix1[1]],data_mod$time)]
if (peak0<peak1) events_mod1$peak[ix1[1]-1] <- events_mod1$peak[ix1[1]]

events_obs1 <- events_obs1[-ix1[1],]
events_mod1 <- events_mod1[-ix1[1],]
}

}

return(list(events_obs_match=events_obs1,
            events_mod_match=events_mod1,
            matched_mod=match_mod,
            matched_obs=match_obs))
}
