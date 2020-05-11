BaseFlowSeperation<-function(streamflow, bf_method, k=0.925, c=quantile(streamflow,probs = 0.25, na.rm=T), filter_parameter=0.925, passes=1){
  if(bf_method=='Duncan'){
    mn<-streamflow
    qflow<-mn
    begin<-1
    end<-length(mn)
    mn_1<-vector(length = end)
    mn[end]<-if(mn[end]<quantile(mn,0.25))mn[end] else mean(mn)/1.5

    for (i in seq((end-1),(begin),by=-1)) {
      mn_1[i] <- ((mn[i+1]-c)/k) + c
      if(mn_1[i]>mn[i]){
        mn_1[i]<-mn[i]
      }

    }
    streamflow<-mn_1
    suppressWarnings(Ends<-c(1,length(streamflow))*rep(1,(passes+1))) # Start and end values for the filter function
    suppressWarnings(AddToStart<-c(1,-1)*rep(1,passes))
    btP<-streamflow##Previous pass's baseflow approximation
    qft<-vector(length=length(streamflow))
    bt<-vector(length=length(streamflow))
    bt[1]<-if(streamflow[1]<quantile(streamflow,0.25)) streamflow[1] else mean(streamflow)/1.5
    ##Guess baseflow value in first time step.
    for(j in 1:passes){
      for (i in (Ends[j]+AddToStart[j]):Ends[j+1]){
        if ((filter_parameter*bt[i-AddToStart[j]]+((1-filter_parameter)/2)*(btP[i]+btP[i-AddToStart[j]]))>btP[i]){
          bt[i]<-btP[i]
        } else bt[i]<-filter_parameter*bt[i-AddToStart[j]]+((1-filter_parameter)/2)*(btP[i]+btP[i-AddToStart[j]])
        qft[i]<-streamflow[i]-bt[i]
      }
      if (j<passes){
        btP<-bt
        bt[Ends[j+1]]<-if(streamflow[Ends[j+1]]<mean(btP))streamflow[Ends[j+1]]/1.2 else mean(btP)
        ##Refines the approximation of end values after the first pass
      }
    }
    f <- data.frame(qflow,bt,qft,mn_1)
    return(f)
  }else if(bf_method=='Lyne-Nathan'){
    suppressWarnings(Ends<-c(1,length(streamflow))*rep(1,(passes+1))) # Start and end values for the filter function
    suppressWarnings(AddToStart<-c(1,-1)*rep(1,passes))
    btP<-streamflow##Previous pass's baseflow approximation
    qft<-vector(length=length(streamflow))
    bt<-vector(length=length(streamflow))
    bt[1]<-if(streamflow[1]<quantile(streamflow,0.25)) streamflow[1] else mean(streamflow)/1.5
    ##Guess baseflow value in first time step.
    for(j in 1:passes){
      for (i in (Ends[j]+AddToStart[j]):Ends[j+1]){
        if ((filter_parameter*bt[i-AddToStart[j]]+((1-filter_parameter)/2)*(btP[i]+btP[i-AddToStart[j]]))>btP[i]){
          bt[i]<-btP[i]
        } else bt[i]<-filter_parameter*bt[i-AddToStart[j]]+((1-filter_parameter)/2)*(btP[i]+btP[i-AddToStart[j]])
        qft[i]<-streamflow[i]-bt[i]
      }
      if (j<passes){
        btP<-bt
        bt[Ends[j+1]]<-if(streamflow[Ends[j+1]]<mean(btP))streamflow[Ends[j+1]]/1.2 else mean(btP)
        ##Refines the approximation of end values after the first pass
      }
    }
    qflow<-streamflow
    f <- data.frame(qflow,bt,qft)
    return(f)
  }else if (bf_method=='Eckhardt'){
    lg <- function(x)c(NA, x[1:(length(x)-1)])
    a<-streamflow
    cc<-data.frame(a)
    cc$b<-lg(a)
    alpha<-lm(cc$a~0+cc$b)$coefficients

    mn<-streamflow
    mn_1<-vector(length = length(mn))
    mn_1[length(mn)]<-if(streamflow[length(mn)]<quantile(streamflow,0.25)) streamflow[length(mn)] else mean(streamflow)/1.5

    for(i in length(mn):2){
      if(mn_1[i]>mn[i]){
        mn_1[i]<-mn[i]
      }
      mn_1[i-1]<-mn_1[i]/alpha
    }
    qflow<-streamflow
    bt<-mn_1
    qft<-qflow-bt
    f <- data.frame(qflow,bt,qft)
    return(f)
  }
}


