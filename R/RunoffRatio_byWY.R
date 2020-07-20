RunoffRatio_byWY<-function(dfRR,WS_area){
  # dfRR: 1st column - streamflow; 2nd column - precip; 
  # 3rd column: Date in POSIXct format

  require(data.table)

  dt1 <- as.data.table(dfRR)
  names(dt1) <- c("Q","PCP","Date")
  dt1[, Qmm:=Q*3600/WS_area*1000]
  dt1[, year:=as.integer(format(Date,"%Y"))]
  dt1[, month:=as.integer(format(Date,"%m"))]
  dt1[, water_year:=ifelse(month<=9,year,year+1)]
  dt1 <- dt1[,c("Qmm","PCP","water_year"),with=F]
  
  dt1[,runoff_ratio:=sum(Qmm)/sum(PCP),by=water_year]
  dt1 <- dt1[,c("water_year","runoff_ratio"),with=F]
  dt1 <- as.data.frame(dt1[!duplicated(dt1),])
  
  return(dt1)
}

