StreamflowElasticity <- function(dfRR, WS_area){

  # dfRR: 1st column - streamflow; 2nd column - precip;
  # 3rd column: Date in POSIXct format

  require(data.table)

  dt1 <- data.table::as.data.table(dfRR)
  print(is.data.table(dt1))

  names(dt1) <- c("Q","PCP","Date")
  dt1[, Qmm:=Q*3600/WS_area*1000]
  dt1[, year:=as.integer(format(Date,"%Y"))]
  dt1[, month:=as.integer(format(Date,"%m"))]
  dt1[, water_year:=ifelse(month<=9,year,year+1)]
  dt1 <- dt1[,c("Qmm","PCP","water_year"),with=F]
  dt1[,sumQ:=sum(Qmm),by=water_year]
  dt1[,sumP:=sum(PCP),by=water_year]
  dt1 <- dt1[,c("sumQ","sumP","water_year"),with=F]
  dt1 <- dt1[!duplicated(dt1),]

  dQ <- diff(dt1$sumQ)
  dP <- diff(dt1$sumP)
  mQ <- mean(dt1$sumQ)
  mP <- mean(dt1$sumP)

  se <- quantile(dQ/dP * (mQ/mP),probs=0.5) 
  return(se)
}   

