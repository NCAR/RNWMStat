\name{BaseFlowIndex}
\alias{BaseFlowIndex}
\title{Hydrograph Baseflow Index}
\usage{
BaseFlowIndex(streamflow, baseflow, na.rm=T)
}
\description{
Estimates baseflow index from streamflow and baseflow time series.
}
\examples{
data("USGS_Streamflow")

#bf_method = 'Duncan'
df<-BaseFlowSeperation(streamflow, bf_method='Duncan', k=0.93, c=quantile(streamflow,probs = 0.25, na.rm=T), filter_parameter=0.925, passes=1)

#Baseflow Index
BFI<-BaseFlowIndex(df$qflow, df$bt, na.rm=T)
}
\arguments{
\item{streamflow}{Streamflow time series}

\item{baseflow}{Baseflow time series }

\item{na.rm}{Logical TRUE or FALSE}
}
\value{
Baseflow Index (range 0-1).
}


