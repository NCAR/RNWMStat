RSquared <- function (model, observation) {
  df<-data.frame(model,observation)
  df2<-na.omit(df)
  lmr2<-lm(observation ~ model, data=df2)
  return(summary(lmr2)$r.squared)
}

