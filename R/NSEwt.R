NseWt <- function (model, observation, w=0.5, p=1) {
  nse   <- hydroGOF::NSE(model, observation, na.rm=TRUE, FUN=NULL, epsilon="Pushpalatha2012")
  lnnse <- hydroGOF::NSE(model, observation, na.rm=TRUE, FUN=log, epsilon="Pushpalatha2012")
  #Weighted mean
  res <- ((w^p) * (nse^p) + (w^p) * (lnnse^p))^(1/p)
  return(res)
}
