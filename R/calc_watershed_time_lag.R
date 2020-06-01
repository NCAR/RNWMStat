# Function for calculating the watershed lag time via cross correlation 
# analysis using basin precipitation and streamflow time series
# Reference: Yilmaz et al., A process-based diagnostic approach to 
#   model evaluation: Application to the NWS distributed hydrologic model, WRR 2008
#
# Questions & suggestions: Yuqiong.Liu@noaa.gov
#
calcLagTime <- function(pcp,flow,prob=0.8,threshold=5.0,max_shift=5) {

###### arguments
# pcp - mean area precipitation time series (hourly)
# flow - streamflow time series (hourly)
# prob - probability threshold for truncating the hourly streamflow time series
#    to focus on high flow events
# threshold - a norminal threshold value (in units of the flow values) for 
#    truncating the time series; set prob to NA if you would like to this type 
#    of threshold instead of the probability threshold
# max_shift - maximum shift in days for calculating the shifted cross correlation
#    5 days should be long enough for slow responding headwater basins
#    while ignoring events driven by seasonal snow.

# make sure pcp and flow time series have the same length
if (length(pcp) != length(flow))
  stop(paste0("ERROR: precipitaiton and flow time seriese must have the same length!"))

# truncate the data based on the threshold
if (!is.na(prob)) {
ix1 <- flow >= quantile(flow,prob)
} else if (!is.na(threshold)) {
ix1 <- flow >= threshold
} else {
stop(paste0("ERROR: either a probability or a norminal threshold should be provided"))
}
 
flow <- flow[ix1]
pcp <- pcp[ix1]
n1 <- length(flow)

# loop through time shifts and identify the shift with the max correlation
corAll <- NULL
mshift <- max_shift*24
if (mshift > (n1-1)) mshift <- n1-1
for (i1 in 0:mshift) {
  x1 <- pcp[1:(n1-i1)]
  y1 <- flow[(1+i1):n1]
  corAll <- c(corAll,cor(x1,y1))
}

lag1 <- which.max(corAll)-1

return(lag1)
}
