rm(list=ls())

library(RNWMStat)

model <- "schaake"
maxGapFill <- 5 #data gaps (model or obs) shorter than this will be spline filled
snow_basins <- c("06289000","09081600","09492400")

dfBasinFinal <- get(load("data/final_list_basins_attributes_short.Rdata"))
ns1 <- nrow(dfBasinFinal)

statsAll <- data.frame()
for (i in 1:ns1) {

gage1 <- dfBasinFinal$gage[i]
area <- dfBasinFinal[["area_gages2(sq km)"]][i]*1000^2
print(paste0(i," ", gage1))

# load model simulation and obs
file1 <- paste0("data/model/",model,"/",gage1,"_obs_mod.csv")
data1 <- read.csv(file1,header=TRUE)
data1$X <- NULL
data1$Date <- as.POSIXct(as.character(data1$Date),format="%Y-%m-%d %H:00:00",tz="UTC")
data1 <- data.table::as.data.table(data1)

# retrict data to the valid period
data1 <- subset(data1,Date>=as.POSIXct("20081001",format="%Y%m%d") &
   Date<=as.POSIXct("2019093023",format="%Y%m%d%H"))

# fill short data gaps with spline interpolation
dates <- seq(min(data1$Date),max(data1$Date), by="hour")
dates1 <- dates[!dates %in% data1$Date]
if (length(dates1)>0) {
data1 <- rbind(data1,data.table::data.table(Date=dates1,mod=NA,PCP=NA,obs=NA))
data1 <- data1[order(data1$Date),]
data1$mod <- zoo::na.approx(data1$mod,maxgap=maxGapFill,na.rm=FALSE)
data1$obs <- zoo::na.approx(data1$obs,maxgap=maxGapFill,na.rm=FALSE)
}
data1 <- na.omit(data1)

# compute non-event stats
Stats_obj <- RNWM_Stats_nonEvent(as.data.frame(data1),area,bf_method="Eckhardt")

# compute event metrics
gw1 <- snow1 <- slow1 <- FALSE
lag1 <- Stats_obj$Obs_LagTime

# groundwater basins in FL
if (gage1 %in% c("02310947","02479300")) gw1 <- TRUE

# slow basins in IL
if (gage1 %in% c("05584500","05525500")) slow1 <- TRUE

# snow basins
if (gage1 %in% snow_basins) snow1 <- 1

Stats_obj1 <- RNWM_Stats_Event(as.data.frame(data1),area,
  snow=snow1,groundwater=gw1, slow=slow1, lagTime=lag1)
Stats_obj <- c(Stats_obj,Stats_obj1)

system("mkdir -p stats")
#save(Stats_obj,file=paste0("stats/stats_obj_",model,"_",gage1,".Rdata"))

# gather scarlar stats and convert into a data frame (containing all basins)
stats1 <- Stats_obj[sapply(Stats_obj,function(x) length(x)==1)]
stats1$Obs_flowSpline <- stats1$Mod_flowSpline <- NULL
stats1 <- as.data.frame(stats1)
stats1$gage <- gage1

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
}
row.names(statsAll) <- 1:ns1

# if seasonal snow basins, watershed lag time does not make much sense
statsAll$Obs_LagTime[statsAll$gage %in% snow_basins] <- NA
statsAll$Mod_LagTime[statsAll$gage %in% snow_basins] <- NA

save(statsAll,file=paste0("stats/stats_all_",model,".Rdata"))
