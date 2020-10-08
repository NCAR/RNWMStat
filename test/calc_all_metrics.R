rm(list=ls())

library(RNWMStat)
source("RNWM_Stats_Event.R")
source("eventIdentification.R")
source("matchEvents.R")

snow_basins <- c("06289000","09081600","09492400","11264500")
flashy_basins <- c("08103900")
slow_basins <- c("02310947")
maxGapFill <- 5
models <- c("SCHAAKE","VIC","XAJ","DVIC_SP","DVIC_GA","DVIC_PH")

dfBasinFinal <- get(load("data/final_list_basins_attributes_short.Rdata"))
ns1 <- nrow(dfBasinFinal)

statsAll <- data.frame()
for (i in 1:ns1) {

gage1 <- dfBasinFinal$gage[i]
area <- dfBasinFinal[["area_sqkm"]][i]*1000^2
print(paste0(i," ", gage1))

snow1 <- FALSE
if (gage1 %in% snow_basins) snow1 <- TRUE
slow1 <- FALSE
if (gage1 %in% slow_basins) slow1 <- TRUE
prob_peak <- 0.9
if (gage1 %in% flashy_basins) prob_peak <- 0.98

k1 <- 0
for (model in models) {

print(paste0(gage1," ", model))

k1 <- k1+1
if (k1==1) {
# load model simulation and obs
data1 <- get(load(paste0("data/calib_valid/",gage1,"_streamflow_precip.RData")))
data1 <- data.table::as.data.table(data1)
names(data1)[names(data1)=="POSIXct"] <- "Date"
# retrict data to the valid period
data1 <- subset(data1,Date>=as.POSIXct("20081001",format="%Y%m%d") &
   Date<=as.POSIXct("2019093023",format="%Y%m%d%H"))
}

data2 <- data1[,c("Date","OBS",model,"PCP"),with=F]
names(data2) <- c("Date","obs","model","PCP")

# fill short data gaps
dates <- seq(min(data2$Date),max(data2$Date), by="hour")
dates1 <- dates[!dates %in% data2$Date]
if (length(dates1)>0) {
tmp <- data.frame(Date=dates1,obs=NA,model=NA,PCP=NA)
data2 <- rbind(data2,tmp)
data2 <- data2[order(Date),]
data2$mod <- zoo::na.approx(data2$mod,na.rm=FALSE,maxgap=maxGapFill)
data2$obs <- zoo::na.approx(data2$obs,na.rm=FALSE,maxgap=maxGapFill)
data2 <- subset(data2,!is.na(mod) & !is.na(obs))
}

# compute non-event stats
Stats_obj <- RNWM_Stats_nonEvent(as.data.frame(data2),area,bf_method="Eckhardt")

# compute event metrics
lag1 <- Stats_obj$Obs_LagTime
thresh1 <- quantile(data2$obs,prob_peak)
Stats_obj1 <- RNWM_Stats_Event(as.data.frame(data2),area,
  snow=snow1,slow=slow1, threshold=thresh1,lagTime=lag1)
Stats_obj <- c(Stats_obj,Stats_obj1)

system("mkdir -p stats")
save(Stats_obj,file=paste0("stats/stats_obj_",model,"_",gage1,".Rdata"))

# gather scarlar stats and convert into a data frame (containing all basins)
stats1 <- Stats_obj[sapply(Stats_obj,function(x) length(x)==1)]
stats1$Obs_flowSpline <- stats1$Mod_flowSpline <- NULL
stats1 <- as.data.frame(stats1)
stats1$gage <- gage1
stats1$model <- model

# reorder columns
strs <- names(stats1)
strs1 <- strs[! substr(strs,1,4) %in% c("Obs_","Mod_")]
strs2 <- strs[substr(strs,1,4) %in% c("Obs_","Mod_")]
strs2 <- gsub("Obs_","",strs2)
strs2 <- gsub("Mod_","",strs2)
strs2 <- unique(strs2)
strs2 <- paste0(rep(c("Obs_","Mod_"),times=length(strs2)),rep(strs2,each=2))
namesInOrder <- c("gage",strs1[strs1!="gage"],strs2)
stats1 <- stats1[,namesInOrder]

statsAll <- rbind(statsAll, stats1)
}}

# if seasonal snow basins, watershed lag time does not make much sense
statsAll$Obs_LagTime[statsAll$gage %in% snow_basins] <- NA
statsAll$Mod_LagTime[statsAll$gage %in% snow_basins] <- NA

save(statsAll,file=paste0("stats/stats_all_",model,".Rdata"))
