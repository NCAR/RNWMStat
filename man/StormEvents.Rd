\name{StormEvents}
\alias{StormEvents}
\title{Number of Storm Events in a Hydrograph}
\usage{
StormEvents(quickflow,qlimit=0)
}
\description{
Estimates number of storm events based on quickflow time series.
}
\examples{
data("USGS_Streamflow")

#Baseflow estimate
df<-BaseFlowSeperation(streamflow, bf_method='Duncan', k=0.93, c=quantile(streamflow,probs = 0.25, na.rm=T), filter_parameter=0.925, passes=1)

#Storm Events
nev<-StormEvents(df$qft,qlimit = 2)

#Plot
plot(nev$quickflow,type = 'l')
lines(nev$nEvents,col='red')
}
\arguments{
\item{quickflow}{Quickflow time series}

\item{qlimit}{Minimum quickflow threshold. The StormEvents() function only consider events higher than this threshold while counting events }

}
\value{
Retun a dataframe with quickflow and corresponding storm event number.
}


