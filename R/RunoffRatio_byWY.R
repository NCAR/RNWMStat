RunoffRatio_byWY<-function(dfRR,WS_area){
  # dfRR: 1st column - streamflow; 2nd column - precip; 
  # 3rd column: Date in POSIXct format

  dt1 <- data.table::as.data.table(dfRR)
  names(dt1) <- c("Q","PCP","Date")
  data.table::dt1[, Qmm:=Q*3600/WS_area*1000]
  data.table::dt1[, year:=as.integer(format(Date,"%Y"))]
  data.table::dt1[, month:=as.integer(format(Date,"%m"))]
  data.table::dt1[, water_year:=ifelse(month<=9,year,year+1)]
  dt1 <- data.table::dt1[,c("Qmm","PCP","water_year"),with=F]
  
  data.table::dt1[,runoff_ratio:=sum(Qmm)/sum(PCP),by=water_year]
  dt1 <- data.table::dt1[,c("water_year","runoff_ratio"),with=F]
  dt1 <- as.data.frame(dt1[!duplicated(dt1),])
  
  return(dt1)
}

