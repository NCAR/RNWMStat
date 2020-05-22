\name{NSE}
\alias{NSE}
\title{NSE (Nash-Sutcliffe efficiency)}
\usage{
NSE(model, observation, nullModel=mean(observation, na.rm=na.rm), na.rm=TRUE)
}
\description{
Nash-Sutcliffe efficiency- between model and observation
}
\examples{
data("USGS_Streamflow")

NSE(model, observation, nullModel=mean(observation, na.rm=na.rm), na.rm=TRUE)

}
\arguments{
\item{model}{model streamflow time series}

\item{observation}{observation streamflow time series }

\item{nullModel}{mean(observation, na.rm=na.rm) }

}
\value{
Retun Nash-Sutcliffe efficiency between model and observation.
}


