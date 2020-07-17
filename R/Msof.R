# multi-scale objective function (MSOF)
# There is not limit on the number of scales to be considered
# The value of scales is defined as the number of time steps of data.
# For example, given daily data, scales=c(1,10,30) means calculate
# MSOF based on daily, 10-day, and 30-day time scales; if it is hourly
# data, scales=c(1,24,72,240,720) means calculate MSOF based on hourly,
# daily, 3-day, 10-day, and 30-day time scales

Msof <- function(m,o, scales=c(1,24,72,120))  {

   if (sum(scales<1)>0) stop("Scales (number of time steps) must not be less than 1!")

   ns <- length(scales)
   n1 <- length(m)
   n2 <- n1 - n1 %% scales
   sum0 <- 0
   for (i in 1:ns) {
     m1 <- m[1:n2[i]]; o1 <- o[1:n2[i]]
     if (scales[i]==1) {
        m2<-m1; o2<-o1
     } else {
        # compute model and observation at the prescribed scales
        m2 <- colMeans(matrix(m1,nrow=scales[i]),na.rm=TRUE)
        o2 <- colMeans(matrix(o1,nrow=scales[i]),na.rm=TRUE)
     }
     # remove missing values in the averaged time series
     # before computing objevtive function
     idx <- !is.na(m2) & !is.na(o2)
     m2 <- m2[idx]; o2 <- o2[idx]

     sum0 <- sum0 + sum((m2-o2)^2)*var(o, na.rm=TRUE)/var(o2)
   }
   obj <- sqrt(sum0)

}

