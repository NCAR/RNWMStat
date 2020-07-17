hyperResMultiObj <- function(m, o, nullModel=mean(o, na.rm=na.rm), na.rm=TRUE) {
   # Input arguments:
   # m - Vector of modeled streamflow (cms)
   # o - Vector of observed streamflow (cms)

   # Establish weights for various metrics
   w0 = 0.4
   w1 = 0.2
   w2 = 0.4

   # First calculate the Normalized Nash Sutcliffe Efficiency (nnse)
   err1 <- sum((m - o)^2, na.rm=na.rm)
   err2 <- sum((o - nullModel)^2, na.rm=na.rm)
   nse <- 1 - (err1/err2)
   NNSE <- 1/(2-nse)

   # Second, calculate the peak discharge error
   Pe <- (max(m,na.rm=TRUE) - max(o,na.rm=TRUE)) / max(o,na.rm=TRUE)

   # Third, calculate the volume error.
   # We are assuming streamflow is cubic meters / second....
   # Also assuming the passed in timestep of observations is in seconds.
   Ve <- sum((m - o), na.rm=TRUE) / sum(o, na.rm=TRUE)

   # Apply a weighting to calulate a final metric.
   objMetric <- w0*(1.0 - NNSE) + w1*abs(Pe) + w2*abs(Ve)
   objMetric
}
