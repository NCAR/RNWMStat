LagTime<-function(pcp,qflow,lag.max=7){

   ccfv<-ccf(pcp,qflow,lag.max = 7,plot = FALSE)
   lagv<-ccfv$lag[,1,1]
   corv<-ccfv$acf[,1,1]
   lTdf<-data.frame(lagv,corv)
   lTdf$lagv[lTdf$lagv<0]<-NA
   lTdf<-na.omit(lTdf)
   return(which(lTdf$corv==max(lTdf$corv,na.rm=T))-1)
}
