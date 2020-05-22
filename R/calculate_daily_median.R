# Calculates the daily median value
#
#--------- FUNCTION SUMMARY -------------
# 1. CalcDailyMed
#    This function calculates the daily median value by
#    aggregating the variable of interest to the daily time step
#    by taking the mean of the values in a day and then calculates
#    the median.
#    REQUIRED INPUTS
#       1. Dataframe with dates and variable of interest (dataframe)
#       2. The column name for the variable of interest (string)
#    OPTIONAL INPUTS
#       1. The column name of the date column (string)
#    DEPENDENCIES
#       None
#----------------------------------------

CalcDailyMed <- function(dataFrame, colName_VariableOfInterest, colName_Date = 'date') {
  names(dataFrame)[names(dataFrame) == colName_VariableOfInterest] <- 'value'
  names(dataFrame)[names(dataFrame) == colName_Date] <- 'date'
  
  # Add a column that includes the month-day-year date format which will be used to
  # aggregate the data over.
  
  dataFrame$date <- format(as.POSIXct(dataFrame$date), format = '%m-%d-%Y')

  
  # Calculate the mean daily values
  
  dataFrame_meanDailyValues <- aggregate(value ~ date, dataFrame, mean)
  
  # Calculate the median of the daily values
  medianDailyValue <- median(dataFrame_meanDailyValues$value)
  
  medianDailyValue

}
