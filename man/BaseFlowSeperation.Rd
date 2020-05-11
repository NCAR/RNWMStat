\name{BaseFlowSeperation}
\alias{BaseFlowSeperation}
\title{Hydrograph Baseflow Seperation}
\usage{
BaseFlowSeperation(streamflow, bf_method, k=0.93, c=quantile(streamflow,probs = 0.25, na.rm=T), filter_parameter=0.925, passes=1))
}
\description{
Estimates baseflow and quickflow from daily streamflow time series using Duncan (2019) method and Lyne and Hollick (1979)-Nathan and McMahon (1990) method.
}
\examples{
data("USGS_Streamflow")

#bf_method = 'Duncan'
Baseflow<-BaseFlowSeperation(streamflow, bf_method='Duncan', k=0.93, c=quantile(streamflow,probs = 0.25, na.rm=T), filter_parameter=0.925, passes=1)

# Plotting
plot(Baseflow$qflow,type = 'l')# Streamflow
lines(Baseflow$bt,col='red') # Baseflow
lines(Baseflow$qft,col='blue') # Quickflow

#bf_method = 'Lyne-Nathan'
Baseflow<-BaseFlowSeperation(streamflow, bf_method='Lyne-Nathan', k=0.93, c=quantile(streamflow,probs = 0.25, na.rm=T), filter_parameter=0.925, passes=1)

# Plotting
plot(Baseflow$qflow,type = 'l')# Streamflow
lines(Baseflow$bt,col='red') # Baseflow
lines(Baseflow$qft,col='blue') # Quickflow
}
\arguments{
\item{streamflow}{Daily streamflow time series}

\item{bf_method}{Baseflow seperation methods (e.g., "Duncan" or "Lyne-Nathan") }

\item{k}{Recession constant. Only for bf_method='Duncan'. Default value is 0.925 based on Nathan and McMahon (1990) }

\item{c}{c is a constant. Only for bf_method='Duncan'. With c positive it simulates a constant flow (such as melt from permanent ice
fields) added to the exponentially receding flow. With c negative it allows modelling of ephemeral streams, by assuming that an exponential decay component greater than zero is still notionally present, but includes some unmeasured hyporheic flow. Default value is 25th quantile of streamflow. }

\item{filter_parameter}{Digital filter parameter for smoothing master recession curve. Default value is 0.925 }

\item{passes}{Number of filter passes. Default value is 1 based on Duncan (2019)}
}
\value{
A dataframe consisting streamflow(qflow), baseflow(bt), quickflow(qft), and master recession (mn_1, available only for bf_method='Duncan' )
}
