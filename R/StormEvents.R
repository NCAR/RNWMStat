StormEvents<-function(quickflow,qlimit=0){
  nEvents<-vector(length = length(quickflow))
  nEvents<-0
  countn <-0
  for (i in 1:length(quickflow)) {
    if(i==1){
      if( (quickflow[i]<=qlimit)&(quickflow[i]>0) ){
        countn <- countn+0
        nEvents[i]<-0
      }else if(quickflow[i]>qlimit){
        countn<-countn+1
        nEvents[i]<-countn
      }
    }else{
      if( (quickflow[i]<=qlimit)&(quickflow[i]>0) ){
        countn <- countn+0
        nEvents[i]<-0
      }else if( (quickflow[i]>qlimit)&(quickflow[i-1]<=qlimit) ){
        countn<-countn+1
        nEvents[i]<-countn
      }else if( (quickflow[i]>qlimit)&(quickflow[i-1]>qlimit) ){
        countn<-countn
        nEvents[i]<-countn
      }else{
        nEvents[i]<-0
      }
    }
  }
  df<-data.frame(quickflow,nEvents)
  return(df)
}
